#' Expression_rate
#' ration between the number of expressed genes and the total number of reads
#' @param genes_by_cells real data
#' @param ngenes_min min number of genes in each cell
#' @return expRate_vect vector with expression rate for each cell
#' @export
 
exp_rate <- function(genes_by_cells, ngenes_min = 5){
  
	cat("Calculating expression rate\n")
  ans <- apply(genes_by_cells, 2, function(x) ifelse(sum(x[x > ngenes_min])>0, sum(x>ngenes_min) / sum(x[x>ngenes_min]), 0))
  names(ans) <- colnames(genes_by_cells)
  
  return(ans)
}