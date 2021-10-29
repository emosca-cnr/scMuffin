#' Calculate the overlap coefficient between two sets
#' @param x vector of items
#' @param y vector of items
#' @export
overlap_coefficient <- function(x, y){
	
	ans <- sum(x %in% y)
	ans <- ans / min(length(x), length(y))
	return(ans)
	
}