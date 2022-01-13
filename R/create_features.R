#' Create feature list
#' @param x data.frame of 1 or more cell values with cell id as row names
#' @export

create_features <- function(x=NULL){
	
  if(!is.data.frame(x)){
    stop("x must be a data.frame")
  }
  
	ans <- list(
	  df=x,
	  type=sapply(x, class)
	)
	
	#non numeric features are converted to factors
	for(i in 1:ncol(ans$df)){
	  if(ans$type[i] == "factor"){
	    ans$df[, i] <- as.factor(ans$df[, i])
	  }
	}
	
	return(ans)
	
}