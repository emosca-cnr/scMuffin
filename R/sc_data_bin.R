#' sc_data_bin
#' 
#' Split the input data matrix by rows into nbins
#' 
#' @param genes_by_cells genes-by-cells matrix
#' @param nbins number of bins
#' @param use.log use logarithm
#' @return vector with bin id for each row
#' @description Split the input data matrix by rows into nbins
#' @export


sc_data_bin <- function(genes_by_cells, nbins=25, use.log=TRUE){


  ans <- rowSums(genes_by_cells)
  if(use.log){
    if(any(ans==0)){
      halfmin <- min(ans[ans>0]) / 2
      ans[ans==0] <- halfmin
    }
    ans <- cut(log2(ans), nbins, labels = FALSE)
  }else{
    ans <- cut(ans, nbins, labels = FALSE)
  }

  return(ans)

}
