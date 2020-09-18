#' gene_set_score
#' @param gene_set_data real data
#' @param control_set_data control data
#' @param nmark_min min number of markers
#' @param ncells_min min number of cells

gene_set_score <- function(gene_set_data, control_set_data, nmark_min=5, ncells_min=NULL){

  #old version
  #ans <- data.frame(case=colMeans(gene_set_data), control=colMeans(control_set_data), case.N=nrow(gene_set_data), case.AV=colSums(gene_set_data!=0), case.NA=colSums(gene_set_data==0), control.N=nrow(control_set_data), control.AV=colSums(control_set_data != 0), control.NA=colSums(control_set_data==0), stringsAsFactors = F)
  #ans$score <- ans$case - ans$control

  #new version
  ans <- list(gs=data.frame(case=colMeans(gene_set_data), case.N=nrow(gene_set_data), case.AV=colSums(gene_set_data!=0), case.NA=colSums(gene_set_data==0), stringsAsFactors = F)) #real

  #z-scores....
  temp <- gene_set_data
  for(i in 1:nrow(temp)){
    temp_idx <- temp[i, ] != 0
    temp[i, temp_idx] <- scale(temp[i, temp_idx])
  }

  ans[[1]] <- cbind(ans[[1]], case.z=colMeans(temp))

  ans <- c(ans, lapply(control_set_data, function(x) data.frame(control=colMeans(x), control.N=nrow(x), control.AV=colSums(x != 0), control.NA=colSums(x==0), stringsAsFactors = F))) #controls

  ans <- lapply(ans, function(x) cbind(x, nmark_min=x[,3]>=nmark_min)) #only cells in which the gene set is available

  #cells with enough genes
  vect_selected <- do.call(cbind, t(lapply(ans, function(x) x$nmark_min))) #the first column contains real data

  #null_ok <- rowSums(vect_selected[, -1])  > (ncol(vect_selected[, -1])) #old
  gene_set_ok_in_nulls <- rowSums(vect_selected[, -1]) #new
  cat("gene_set_ok_in_nulls...\n")
  print(table(gene_set_ok_in_nulls[gene_set_ok_in_nulls>0]))

  null_ok <- gene_set_ok_in_nulls > (ncol(vect_selected)-1)*0.9 #new

  if(sum(vect_selected[, 1]) >= ncells_min){

    vect_selected[, 1] <- vect_selected[, 1] & null_ok

    if(sum(vect_selected[, 1]) >= ncells_min){

      #score matrix
      score <- do.call(cbind, t(lapply(ans[-1], function(x) ans[[1]]$case - x$control)))
      score <- score * sign(vect_selected[, -1])
      score <- rowMeans(score)

      #average control value for each cell
      avg_control <- do.call(cbind, t(lapply(ans[-1], function(x) x$control)))
      avg_control <- avg_control * sign(vect_selected[, -1])
      avg_control <- rowMeans(avg_control)

      res <- data.frame(ans[[1]], avg_control=avg_control, avg_delta_score=score, stringsAsFactors = F) #difference between real gs score (1), and each of the null values [2, k]

      #res <- res[res$nmark_min, ] #eliminate cells without a sufficient number of genes

    }else{
      warning("cannot create null models for at least ", ncells_min, " cells\n")
      res <- data.frame(ans[[1]], avg_control=NA, avg_delta_score=NA, stringsAsFactors = F)
    }

    res <- res[res$nmark_min, ] #eliminate cells without a sufficient number of genes

  }else{

    warning("cannot calcolate gene set score for at least ", ncells_min, " cells\n")
    res <- NA

  }

  ans <- list(gs_score=ans, res=res, cells_ok=vect_selected)

  return(ans)

}
