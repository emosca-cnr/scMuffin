#' Expression_rate
#' @param genes_by_cells real data
#' @param ngenes_min min number of genes in each cell
#' @return expRate_vect vector with expression rate for each cell
 
exp_Rate <- function(genes_by_cells, ngenes_min = 5){
  
  expRate_vect <- apply(genes_by_cells,2, function(x) sum(x>5)/sum(x[x>5]))
  names(expRate_vect) <- colnames(genes_by_cells)
  
  return(expRate_vect)
}