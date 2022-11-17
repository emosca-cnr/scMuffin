#' Transcriptional Complexity
#' @param scMuffinList scMuffinList object
#' @param min_counts minimum number of counts
#' @param min_cells minimum number of cells in which a genes must have at least min_counts counts
#' @param min_genes minimum number of genes that a cell must have with at least min_counts counts
#' @return scMuffinList with transcr_compl element
#' @export

transcr_compl <- function(scMuffinList = NULL, min_counts = 5, min_cells=10, min_genes=500){

  if(length(scMuffinList$counts)==0){
	  stop("scMuffinList does not contain counts\n")
  }
  counts <- scMuffinList$counts

  cat("Cleaning\n")
  counts[counts < min_counts] <- 0
  counts[Matrix::rowSums(sign(counts)) < min_cells, ] <- 0 #only genes in at least 10 cells
  counts[, Matrix::colSums(sign(counts)) < min_genes] <- 0 #only cells in with least 100 genes
 
  N <- Matrix::colSums(counts)

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
    
  scMuffinList$transcr_compl <- list(summary=data.frame(tot_counts=as.numeric(Tot_transcripts), n_genes=as.numeric(ngenes), C=as.numeric(tr_compl), H=as.numeric(ans), row.names=names(ans)), full=list())

  return(scMuffinList)

}

