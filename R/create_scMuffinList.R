#' Create scMuffinList
#' @description Create scMuffinList from counts and/or normalized data. 
#' @param counts genes-by-cell matrix or data.frame of counts.
#' @param normalized genes-by-cell matrix or data.frame with normalized expression.
#' @return scMuffinList with counts and normalized elements as dgCMatrix.
#' @importFrom Matrix Matrix

create_scMuffinList <- function(counts=NULL, normalized=NULL){
  
  scMuffinList <- list()
  
  if(!is.null(counts) | !is(counts, "dgCMatrix")){
    scMuffinList$counts <- Matrix::Matrix(counts, sparse = TRUE)
    cat("counts created\n")
  }

  if(!is.null(normalized) | !is(normalized, "dgCMatrix")){
    scMuffinList$normalized <- Matrix::Matrix(normalized, sparse = TRUE)
    cat("normalized created\n")
  }
  
  return(scMuffinList)
}