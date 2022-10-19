#' Create feature list
#' @param x data.frame of 1 or more cell values with cell id as row names
#' @description Create feature list
#' @export

create_features_obj <- function(x=NULL){
	
  if(!is.data.frame(x)){
    stop("x must be a data.frame")
  }
  
	ans <- list(df=x)
	
	#non numeric features are converted to factors
	for(i in 1:ncol(ans$df)){
	  if(!is.numeric(ans$df[, i])){
	    ans$df[, i] <- as.factor(ans$df[, i])
	  }
	}
	
	ans$type <- sapply(ans$df, class)
	
	return(ans)
	
}