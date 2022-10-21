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
  if(na.rm){
    n <- Matrix::rowSums(genes_by_cells>0)
    ans <- ans/n
  }
  #if(use.log){
    # if(any(ans==0)){
    #   halfmin <- min(ans[ans>0]) / 2
    #   ans[ans==0] <- halfmin
    # }
  #  ans <- cut(log1p(ans), nbins, labels = FALSE)
  #}else{
  #  ans <- cut(ans, nbins, labels = FALSE)
  #}

  ans <- ggplot2::cut_number(ans, n = nbins, labels=FALSE)
  
  return(ans)

}
