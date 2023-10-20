#' Split the genes into bins by their expression
#' @param genes_by_cells genes-by-cells matrix
#' @param nbins number of bins
#' @param na.rm if TRUE, cells with null (0) values are not considered in the calculation of the gene expression level.
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowSums
#' @return vector with bin identifier for each gene
#' @description Split the input data matrix by rows into nbins
#' @export


sc_data_bin <- function(genes_by_cells=NULL, nbins=25, na.rm=FALSE){


  ans <- rowSums(genes_by_cells)
  
  #if excluding NA every gene will have a different number of cells; therefore mean is used to define the expression bins
  if(na.rm){
    n <- rowSums(genes_by_cells>0)
    ans <- ans/n
  }
 
  ans <- cut_number(ans, n = nbins, labels=FALSE)
  
  return(ans)

}
