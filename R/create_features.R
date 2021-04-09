#' Create feature list
#'
#' @export

create_features <- function(cell_id=NULL, values=NULL){
	
	ans <- list(df=data.frame(values, stringsAsFactors = F, row.names = cell_id))

	ans$type <- unlist(lapply(ans$df, class))
	
	return(ans)
	
}