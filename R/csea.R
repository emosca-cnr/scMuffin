#' Gene Set Enrichment Analysis
#' @param rl numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory
#' @param gsl named list of gene sets
#' @param k number of permutations
#' @param ord.mode ordering mode: -1 -> descending; 1 ascending; must be of length equal to `ncol(rl)`
#' @param mc_cores_path number of cores to use for parallel calculation of gene set lists; the total number of cpu used will be mc_cores_path * mc_cores_perm
#' @param mc_cores_perm number of cores to use for parallel calculation of ranked list permutations; the total number of cpu used will be mc_cores_path * mc_cores_perm
#' @param min.k minimum number of valid permutations to support empirical nulls
#' @param min.size minimum number of cells with a not null value
#' @import parallel
#' @importFrom stats p.adjust
#' @importFrom qvalue qvalue
#' @return list with two data.frames, gs_table and leading_edge. 
#` \enumerate{
#`  \item {gs_table} a data.frame with:
#`  \itemize{
#`   \item es, enrichment score; 
#`   \item nes normalized enrichment score;
#`   \item p-value, empirical p-value;
#`   \item adjusted p-value, BH FDR;
#`   \item FDR q-value, empirical FDR. 
#`  }
#`  \item {leading_edge} contains:
#`  \itemize{
#`   \item tags, leading edge size;
#`   \item tags_perc, leading edge size percent over gene set;
#`   \item list_top, rank of the ES;
#`   \item list_top_perc, rank of the ES percent over full ranked list;
#`   \item lead_edge, gene names of the leading edge.
#`  }
#` }
#' @export


csea <- function(rl=NULL, gsl=NULL, k=100, min.size=100, ord.mode=-1, min.k=10, mc_cores_path=1, mc_cores_perm=1){

  #cheks
  if(!is.matrix(rl) | !is.numeric(rl)){
    stop("rl must be a numeric matrix")
  }

  if(length(ord.mode) != ncol(rl)){
    stop("length(ord.mode) must be equal to ncol(rl)")
  }

  min.k <- min(min.k, k)
  
  #create the list of ranked vectors
  rll <- vector('list', ncol(rl))
  names(rll) <- colnames(rl)
  for(i in 1:length(rll)){
    rll[[i]] <- sort(array(rl[, i], dimnames = list(rownames(rl))), decreasing = ord.mode[i]==-1)
  }
  
  ##keep just elements that are !=0
  rll <- lapply(rll, function(x) x[x!=0])
  rll <- rll[lengths(rll)>min.size]
  cat("Ranked list that passed the checks:", names(rll), "\n")
  
  #real es
  cat("ES of input data...\n")
  
  gsl_size <- lengths(gsl)
  
  real_es_data <- lapply(gsl, function(x) lapply(rll, function(y) es(which(names(y) %in% x), y)))
  real_es <- do.call(rbind, lapply(real_es_data, function(x) unlist(lapply(x, function(y) y$es))))
  
  leading_edge <- vector("list", length(rll))
  names(leading_edge) <- names(rll)
  for(i in 1:length(leading_edge)){
    leading_edge[[i]] <- do.call(rbind, lapply(real_es_data, function(x) x[[i]][, -1]))
  }

  #permutations
  #k permutation of all gene ids
  cat('generating', k, 'permutations\n')
  x_perm <- lapply(1:k, function(x) sample(rownames(rl), nrow(rl)))
  
  #every element of x_perm is a list with a rll-specific permutation, because each element of rll has different size
  x_perm <- vector("list", length = k)
  for(i in 1:length(x_perm)){
    x_perm[[i]] <- lapply(rll, function(x) sample(rownames(x), length(x)))
  }
  
  cat("ES of permutations...\n")
  res <- mclapply(gsl, function(x) do.call(rbind, mclapply(x_perm, function(y) calc_gs_perm(rll, y, x), mc.cores=mc_cores_perm)), mc.cores = mc_cores_path)
 
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

    n_pos_perm <- rowSums(res[[i]]>0)
    n_neg_perm <- rowSums(res[[i]]<0)
    
    p_val <- apply(res[[i]], 1, function(x) ifelse(x[1] >= 0, sum(x >= x[1]) / length(x[x>=0]), sum(x <= x[1]) / length(x[x<=0])))
    
    #p-values for real ES == 0 are set to 1
    #p_val[res[[i]][, 1] == 0] <- 1
    
    idx_na <- which(res[[i]][, 1] == 0)
    if(length(idx_na)>0){
      cat("\tES==0:\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      p_val[idx_na] <- NA ### 
    }
    idx_na <- which(res[[i]][, 1] > 0 & n_pos_perm < min.k)
    if(length(idx_na)>0){
      cat("\tES>0 but less than", min.k, " positive ES in permutations\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      p_val[idx_na] <- NA ### 
    }
    idx_na <- which(res[[i]][, 1] < 0 & n_neg_perm < min.k)
    if(length(idx_na)>0){
      cat("\tES<0 but less than", min.k, " negative ES in permutations\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      p_val[idx_na] <- NA ### 
    }
    
    #normalized ES
    means <- t(apply(res[[i]], 1, function(x) c(mean(x[x>0]), abs(mean(x[x<0]))))) #positive, negative
    means[is.nan(means)] <- 0 #NaN values are caused by the absence of any positive or negative value
    nes <- res[[i]] / means[, 1]
    nes_neg <- res[[i]] / means[, 2]
    nes[res[[i]] < 0] <- nes_neg[res[[i]] < 0]
    nes[is.nan(nes)] <- NA #NaN values are caused by 0/0
    rm(means, nes_neg)
    
    #if there are not at least min.k the NES is unrelieable
    nes[nes>0 & n_pos_perm < min.k] <- NA
    nes[nes<0 & n_neg_perm < min.k] <- NA

    #calculate FDR
    all_nes <- as.numeric(nes)
    n_nes_pos <- sum(nes>0, na.rm = T)
    n_nes_neg <- sum(nes<0, na.rm = T)
    n_real_nes_pos <- sum(nes[,1] > 0, na.rm = T)
    n_real_nes_neg <- sum(nes[,1] < 0, na.rm = T)
 
    #FDR: NES* > 0: fdrq = #(all positive NESp >= NES*) / #(all positive NESp) / [ #(all NES* >= NES*) / (all positive NES*) ]
    #FDR: NES* < 0: fdrq = #(all negative NESp <= NES*) / #(all negative NESp) / [ #(all NES* <= NES*) / (all negative NES*) ]
    
    fdrq <- sapply(nes[, 1], function(x) ifelse(x>0, sum(all_nes >= x, na.rm = T) / n_nes_pos, sum(all_nes <= x, na.rm = T) / n_nes_neg))
    fdrq <- fdrq / sapply(nes[, 1], function(x) ifelse(x>0, sum(nes[, 1] >= x, na.rm = T) / n_real_nes_pos, sum(nes[, 1] <= x, na.rm = T) / n_real_nes_neg))
  
    #q of nes that are equal to 0 are set to 1
    #fdrq[nes[, 1] == 0] <- 1

    rm(all_nes)
    fdrq[fdrq>1] <- 1
    idx_replace <- which(fdrq<p_val)
    fdrq[idx_replace] <- p_val[idx_replace]
    
    #### FROM DOSE
    qobj <- tryCatch(qvalue(p_val, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
    if (inherits(qobj, "qvalue")) {
      qvalues <- qobj$qvalues
    } else {
      qvalues <- NA
    }
    
    nperm <- n_pos_perm
    nperm[which(res[[i]][, 1]<0)] <- n_neg_perm[which(res[[i]][, 1]<0)]
    nperm[which(res[[i]][, 1] == 0)] <- NA
    
    #out table
    #out[[i]] <- data.frame(id=rownames(res[[i]]), es=res[[i]][, 1], p_val=p_val, adj_p_val=p.adjust(p_val, method='fdr'), n_pos_perm=n_pos_perm, n_neg_perm=n_neg_perm, nes=nes[, 1], FDRq=fdrq, stringsAsFactors=FALSE)
    
    #out table
    out[[i]] <- data.frame(id=rownames(res[[i]]), size=gsl_size[match(rownames(res[[i]]), names(gsl_size))], es=res[[i]][, 1], nes=nes[, 1], nperm=nperm, p_val=p_val, adj_p_val=p.adjust(p_val, method='fdr'), q_val=qvalues, FDRq=fdrq, stringsAsFactors=FALSE)
    
    out[[i]] <- merge(out[[i]], leading_edge[[i]], by.x=1, by.y=0, sort=F)
    
    out[[i]] <- out[[i]][, -which(colnames(out[[i]]) == "lead_edge_subset")]
    
  }

  print("done")

  #return(list(gs_table=out, leading_edge=leading_edge))
  return(out)

}
