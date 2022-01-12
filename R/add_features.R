#' Adds a feature to the feature list
#' @param features features object
#' @param features_to_add data.frame of 1 or more cell values with cell id as row names
#' @export

add_features <- function(features=NULL, features_to_add=NULL){
	
  ans <- create_features(x = features_to_add)
  
  ans$df <- merge(features$df, ans$df, by=0, all=T, sort=F)
  rownames(ans$df) <- ans$df[, 1]
  ans$df[, 1] <- NULL
	
  ans$type <- c(features$type, ans$type)
	
	return(ans)
	
}