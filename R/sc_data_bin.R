#' sc_data_bin
#' 
#' Split the input data matrix by rows into nbins
#' 
#' @param genes_by_cells genes-by-cells matrix
#' @param nbins number of bins
#' @param na.rm whether to remove NA or not
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowSums
#' @return vector with bin id for each row
#' @description Split the input data matrix by rows into nbins
#' @export


sc_data_bin <- function(genes_by_cells, nbins=25, na.rm=FALSE){


  ans <- Matrix::rowSums(genes_by_cells)
  
  #if excluding NA every gene will have a different number of cells; therefore mean is used to define the expression bins
  if(na.rm){
    n <- Matrix::rowSums(genes_by_cells>0)
    ans <- ans/n
  }
 
  ans <- ggplot2::cut_number(ans, n = nbins, labels=FALSE)
  
  return(ans)

}
