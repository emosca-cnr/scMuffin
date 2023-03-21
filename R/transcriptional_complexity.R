#' Transcriptional Complexity and Entropy
#' @param scMuffinList scMuffinList object
#' @param min_counts minimum number of counts
#' @param min_cells minimum number of cells in which a gene must have at least min_counts counts
#' @param min_genes minimum number of genes that a cell must express with at least min_counts counts
#' @return scMuffinList with transcr_compl element a list with summary (full is empty):
#' \itemize{
#'   \item{tot_counts, total number of transcripts;}
#'   \item{n_genes, total number of expressed genes;}
#'   \item{C, transcriptional complexity;}
#'   \item{H, transcriptional entropy;}
#'}
#' @export

transcr_compl <- function(scMuffinList = NULL, min_counts = 5, min_cells=10, min_genes=500){
  
  if(length(scMuffinList$counts)==0){
    stop("scMuffinList does not contain counts\n")
  }
  counts <- scMuffinList$counts
  
  cat("Cleaning\n")
  idx <- counts < min_counts
  if(any(idx)){
    counts[idx] <- 0
  }
  idx <- Matrix::rowSums(sign(counts)) < min_cells
  if(any(idx)){
    counts[idx, ] <- 0 #only genes in at least 10 cells
  }
  idx <- Matrix::colSums(sign(counts)) < min_genes
  if(any(idx)){
    counts[, idx] <- 0 #only cells in with least 100 genes
  }
  
  N <- Matrix::colSums(counts)
  if(sum(N>0) < min_cells){
    stop("not enough cells after filtering")
  }
  
  #Entropy
  cat("Entropy\n")
  ans <- t(apply(counts, 1, function(i_row) i_row / N))
  ans[is.infinite(ans) | is.nan(ans)] <- 0
  ans <- -apply(ans, 2, function(i_col) sum(i_col*log(i_col+1e-10)))
  
  #Complexity
  cat("Complexity\n")
  ngenes <- apply(counts, 2, function(x) sum(sign(x)))
  Tot_transcripts <- apply(counts, 2, function(x) sum(x))
  tr_compl <- ngenes / Tot_transcripts
  tr_compl[is.infinite(tr_compl) | is.nan(tr_compl)] <- 0
  tr_compl <- tr_compl * max(Tot_transcripts) / max(ngenes)
  
  #Residual of tot_counts ~ n_genes
  cat("Linear model\n")
  xx <- log10(as.numeric(Tot_transcripts))
  yy <- log10(as.numeric(ngenes))
  idx_non_zero <- which(xx>0)
  m <- glm( yy[idx_non_zero] ~ xx[idx_non_zero])
  
  lm <- rep(0, length(xx))
  lm[idx_non_zero] <- m$residuals
  
  scMuffinList$transcr_compl <- list(summary=data.frame(tot_counts=as.numeric(Tot_transcripts), n_genes=as.numeric(ngenes), C=as.numeric(tr_compl), H=as.numeric(ans), LM=lm, row.names=names(ans)), full=list(LM=m))
  
  return(scMuffinList)
  
}

