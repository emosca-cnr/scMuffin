#' Expression_rate
#' number of expressed genes by sum of their counts
#' 
#' Ratio between the number of expressed genes and the total number of reads
#' @param genes_by_cells count matrix
#' @param min_counts min_counts
#' @return expRate_vect vector with expression rate for each cell
#' @export

exp_rate <- function(genes_by_cells, min_counts = 5){
  
	cat("Calculating expression rate\n")
  ans <- apply(genes_by_cells, 2, function(x) ifelse(sum(x > min_counts)>0, sum(x>min_counts) / sum(x[x>min_counts]), 0))
  names(ans) <- colnames(genes_by_cells)
  
  return(ans)
}