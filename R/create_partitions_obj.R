#' Create the clustering object
#' @param x data.frame of 1 or more cell labels with cell id as row names
#' @description Create the clustering object
#' @export

create_partitions_obj <- function(x){
  
  if(!is.data.frame(x)){
    stop("x must be a data.frame")
  }
	
  #clusterings must be factors
  for(i in 1:ncol(x)){
    x[, i] <- as.factor(x[, i])
  }
  
  return(x)
	
}