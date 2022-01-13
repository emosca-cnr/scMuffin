#' Gene Sert Enrichment Analysis
#' @param rl numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory
#' @param gsl named list of gene sets
#' @param k integer, number of permutations
#' @param ord.mode ordering mode: -1 -> descending; 1 ascending; must be of lenght equal to ncol(rl)
#' @param mc_cores_path number of cores to use for parallel calculation of gene set lists; the total number of cpu used will be mc_cores_path x mc_cores_perm
#' @param mc_cores_perm number of cores to use for parallel calculation of ranked list permutations; the total number of cpu used will be mc_cores_path x mc_cores_perm
#' @import parallel
#' @return data.frame with es, nes, p-value, adjusted p-value and FDR q-value
#' @importFrom stats p.adjust

gsea <- function(rl, gsl, k=100, ord.mode=-1, mc_cores_path=1, mc_cores_perm=1){

  #cheks
  if(!is.matrix(rl) | !is.numeric(rl)){
    stop("rl must be a numeric matrix")
  }

  if(length(ord.mode) != ncol(rl)){
    stop("length(ord.mode) must be equal to ncol(rl)")
  }

  #create the list of ranked vectors
  rll <- vector('list', ncol(rl))
  names(rll) <- colnames(rl)
  for(i in 1:length(rll)){
    rll[[i]] <- sort(array(rl[, i], dimnames = list(rownames(rl))), decreasing = ord.mode[i]==-1)
  }

  #permutation of gene ids
  cat('generating', k, 'permutations\n')
  x_perm <- lapply(1:k, function(x) sample(rownames(rl), nrow(rl)))

  #real es
  print("ES...")
  #real_es <- lapply(gsl, function(x) lapply(rll, function(y) es(which(names(y) %in% x), y)))
  #real_es <- do.call(rbind, real_es)

  real_es_data <- lapply(gsl, function(x) lapply(rll, function(y) es(which(names(y) %in% x), y, le=T)))
  real_es <- do.call(rbind, lapply(real_es_data, function(x) unlist(lapply(x, function(y) y$es))))

  leading_edge <- vector("list", length(rll))
  names(leading_edge) <- names(rll)
  for(i in 1:length(leading_edge)){
    leading_edge[[i]] <- do.call(rbind, lapply(real_es_data, function(x) x[[i]]$lea))
  }

  #permutations
  print("calculating permutations")
  if(mc_cores_path==1){
    if(mc_cores_perm == 1){
      #res <- lapply(gsl, function(x) unlist(lapply(x_perm, function(y) calc_gs_perm(rl, y, x))))
      res <- lapply(gsl, function(x) do.call(rbind, lapply(x_perm, function(y) calc_gs_perm(rll, y, x))))
    }else{
      cat(k, "permutations on", mc_cores_perm, "cores\n")
      #res <- lapply(gsl, function(x) unlist(mclapply(x_perm, function(y) calc_gs_perm(rl, y, x), mc.cores=mc_cores_perm)))
      res <- lapply(gsl, function(x) do.call(rbind, parallel::mclapply(x_perm, function(y) calc_gs_perm(rll, y, x), mc.cores=mc_cores_perm)))
    }
  }else{
    cat(length(gsl), "gene sets on", mc_cores_path, "cores\n")
    if(mc_cores_perm == 1){
      #res <- parallel::mclapply(gsl, function(x) unlist(lapply(x_perm, function(y) calc_gs_perm(rl, y, x))), mc.cores = mc_cores_path)
      res <- parallel::mclapply(gsl, function(x) do.call(rbind, lapply(x_perm, function(y) calc_gs_perm(rll, y, x))), mc.cores = mc_cores_path)
    }else{
      cat(k, "permutations on", mc_cores_perm, "cores\n")
      #res <- parallel::mclapply(gsl, function(x) unlist(parallel::mclapply(x_perm, function(y) calc_gs_perm(rl, y, x), mc.cores=mc_cores_perm)), mc.cores = mc_cores_path)
      res <- parallel::mclapply(gsl, function(x) do.call(rbind, parallel::mclapply(x_perm, function(y) calc_gs_perm(rll, y, x), mc.cores=mc_cores_perm)), mc.cores = mc_cores_path)
    }
  }


  #list of pathwa-by-k matrices of permuted es
  #res <- do.call(rbind, res)
  #   if(nrow(res) != length(gsl) | ncol(res) != k){
  #     stop("not all the permutations returned a correct value\n")
  #   }
  #res <- cbind(real_es, res)

  temp <- vector('list', length(rll))
  names(temp) <- colnames(rl)
  for(i in 1:length(rll)){
    temp[[i]] <- cbind(real_es[, i], do.call(rbind, lapply(res, function(x) x[, i])))
  }
  res <- temp
  rm(temp)




  #statistics
  print("calculating statistics...")

  #the first column is the real value, and it is included in the calculation of p
  #if ES* > 0 -> p = # (ESp >= ES*) / (k+1)
  #if ES* < 0 -> p = # (ESp <= ES*) / (k+1)
  out <- res

  for(i in 1:length(rll)){

    p_val <- apply(res[[i]], 1, function(x) ifelse(x[1] >= 0, sum(x >= x[1]) / length(x), sum(x <= x[1]) / length(x)))

    #normalized ES
    means <- t(apply(res[[i]], 1, function(x) c(mean(x[x>0]), abs(mean(x[x<0]))))) #positive, negative
    nes <- res[[i]] / means[, 1]
    nes_neg <- res[[i]] / means[, 2]
    nes[res[[i]] < 0] <- nes_neg[res[[i]] < 0]
    rm(means, nes_neg)

    #calculate FDR
    all_nes <- as.numeric(nes)
    n_nes_pos <- sum(nes>=0)
    n_nes_neg <- sum(nes<=0)
    n_real_nes_pos <- sum(nes[,1] >= 0)
    n_real_nes_neg <- sum(nes[,1] <= 0)
    if((n_nes_pos + n_nes_neg) != (k*length(gsl) + nrow(nes))){
      warning('some nes value is equal to zero')
    }

    #FDR: NES* > 0: fdrq = #(all positive NESp >= NES*) / #(all positive NESp) / [ #(all NES* >= NES*) / (all positive NES*) ]
    #FDR: NES* < 0: fdrq = #(all negative NESp <= NES*) / #(all negative NESp) / [ #(all NES* <= NES*) / (all negative NES*) ]
    fdrq <- sapply(nes[, 1], function(x) ifelse(x>0, sum(all_nes >= x) / n_nes_pos, sum(all_nes <= x) / n_nes_neg))
    fdrq <- fdrq / sapply(nes[, 1], function(x) ifelse(x>0, sum(nes[, 1] >= x) / n_real_nes_pos, sum(nes[, 1] <= x) / n_real_nes_neg))
    rm(all_nes)
    fdrq[fdrq>1] <- 1

    #out table
    out[[i]] <- data.frame(id=rownames(res[[i]]), es=res[[i]][, 1], p_val=p_val, adj_p_val=stats::p.adjust(p_val, method='fdr'), nes=nes[, 1], FDRq=fdrq, stringsAsFactors=FALSE)

  }

  print("done")

  return(list(gs_table=out, leading_edge=leading_edge))

}
