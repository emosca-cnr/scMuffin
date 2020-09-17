#' sc_data_bin
#' @param seurat_data seurat data object
#' @param nbins number of bins
#' @param use.log use logarithm
##' @import libreria_necessaria


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
