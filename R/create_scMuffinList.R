#' Create scMuffinList
#' @description Create scMuffinList from counts and/or normalized data. The gene identifiers must be the same. Note that CNV analysis is currently available only for symbols as gene identifiers.
#' @param counts genes-by-cells matrix or data.frame of counts. Rownames must be gene identifiers. 
#' @param normalized genes-by-cells matrix or data.frame with normalized expression. Rownames must be gene identifiers. 
#' @return scMuffinList with counts and normalized elements as dgCMatrix.
#' @importFrom Matrix Matrix
#' @importFrom methods is
#' @export
#' 
create_scMuffinList <- function(counts=NULL, normalized=NULL){
  
  scMuffinList <- list()
  
  if(!is.null(counts) & !is.null(normalized)){
    if(nrow(counts) != nrow(normalized)){
      stop("counts and normalized matrices have different number of rows.\n")
    }else{
      normalized <- normalized[match(rownames(counts), rownames(normalized)), match(colnames(counts), colnames(normalized))]
      if(nrow(counts) != nrow(normalized)){
        stop("counts and normalized matrices have different number of rows after matching identifiers.\n")
      }
    }
  }
  
  idx_zero <- which(rowSums(counts)==0)
  if(length(idx_zero)>0){
    
    cat("Found genes with all-zero values. These genes may cause issues and will be removed.\n")
    counts <- counts[-idx_zero, ]
    normalized <- normalized[-idx_zero, ]
    
    idx_zero <- which(colSums(counts)==0)
    if(length(idx_zero)>0){
      counts <- counts[, -idx_zero]
      normalized <- normalized[, -idx_zero]
    }
    
    cat("Updated size:\n")
    print(dim(counts))
    
  }
  
  scMuffinList$counts <- Matrix::Matrix(counts, sparse = TRUE)
  cat("counts created\n")
  
  scMuffinList$normalized <- Matrix::Matrix(normalized, sparse = TRUE)
  cat("normalized created\n")
 
   
  
  
  return(scMuffinList)
}