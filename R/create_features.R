#' Create feature list
#' @param cell_id vector of cell id
#' @param values vector of values
#' @export

create_features <- function(cell_id=NULL, values=NULL){
	
	ans <- list(df=data.frame(values, stringsAsFactors = F, row.names = cell_id))

	ans$type <- unlist(lapply(ans$df, class))
	
	return(ans)
	
}