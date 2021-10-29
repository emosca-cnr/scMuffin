#' Create clustering list
#' @param cell_id vector of cell id
#' @param values vector of values
#' @export

create_clusterings <- function(cell_id=NULL, values=NULL){
	
	ans <- data.frame(values, row.names = cell_id)
	
	return(ans)
	
}