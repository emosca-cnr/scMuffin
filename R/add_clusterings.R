#' Adds a clustering to the clustering list
#' @param clusterings clustering list
#' @param cell_id vector of cell identifiers
#' @param values vector or data.frame with clustering labels
#' @export

add_clusterings <- function(clusterings=NULL, cell_id=NULL, values=NULL){
	
	for(i in 1:ncol(values)){
		if(is.integer(values[, i])){
			values[, i] <- as.character(as.numeric(values[,i]))
		}
	}
	temp <- data.frame(values, row.names = cell_id)

	clusterings <- merge(clusterings, temp, by=0, all=T, sort=F)
	rownames(clusterings) <- clusterings[, 1]
	clusterings[, 1] <- NULL
	
	return(clusterings)
	
}