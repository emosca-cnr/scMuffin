#' Create scMuffinList
#' @description Create scMuffinList from counts and/or normalized data. The gene identifiers must be the same. Note that CNV analysis is currently available only for symbols as gene identifiers.
#' @param counts genes-by-cells matrix or data.frame of counts. Rownames must be gene identifiers. 
#' @param normalized genes-by-cells matrix or data.frame with normalized expression. Rownames must be gene identifiers. 
#' @return scMuffinList with counts and normalized elements as dgCMatrix.
#' @importFrom Matrix Matrix
#' @importFrom methods is
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