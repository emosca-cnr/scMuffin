#' Create the clustering object
#' @param x data.frame of 1 or more cell labels with cell id as row names
#' @export

create_clusterings <- function(x){
  
  if(!is.data.frame(x)){
    stop("x must be a data.frame")
  }
	
  #clusterings must be factors
  for(i in 1:ncol(x)){
    x[, i] <- as.factor(x[, i])
  }
  
  return(x)
	
}