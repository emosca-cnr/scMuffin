#' Adds a feature to the feature list
#' @param clusterings feature list
#' @param cell_id vector of cell identifiers
#' @param values vector or data.frame with feature values
#' @export

add_features <- function(features=NULL, cell_id=NULL, values=NULL){
	
	temp <- data.frame(values, stringsAsFactors = F, row.names = cell_id)

	features$df <- merge(features$df, temp, by=0, all=T, sort=F)
	rownames(features$df) <- features$df[, 1]
	features$df[, 1] <- NULL
	
	features$type <- unlist(lapply(features$df[-1], class))
	
	return(features)
	
}