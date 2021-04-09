#' Create feature list
#'
#' @export

create_clusterings <- function(cell_id=NULL, values=NULL){
	
	ans <- data.frame(values, row.names = cell_id)
	
	return(ans)
	
}