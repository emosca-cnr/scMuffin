#' Adds a feature to the feature list
#' @param features_obj features object
#' @param to_add data.frame of 1 or more cell values with cell id as row names
#' @export

add_features <- function(features_obj=NULL, to_add=NULL){
	
  ans <- create_features_obj(x = to_add)
  
  ans$df <- merge(features_obj$df, ans$df, by=0, all=T, sort=F)
  rownames(ans$df) <- ans$df[, 1]
  ans$df[, 1] <- NULL
	
  ans$type <- c(features_obj$type, ans$type)
	
	return(ans)
	
}