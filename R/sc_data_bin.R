#' sc_data_bin
#' 
#' Split the input data matrix by rows into nbins
#' 
#' @param seurat_data seurat data object
#' @param nbins number of bins
#' @param use.log use logarithm
#' @value vector with bin id for each row


sc_data_bin <- function(seurat_data, nbins=25, use.log=TRUE){


  ans <- rowSums(seurat_data)
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
