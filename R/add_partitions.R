#' Adds a clustering to the clustering list
#' @param partitions_obj clusterings object
#' @param to_add data.frame of 1 or more cell labels with cell id as row names
#' @description Adds a clustering to the clustering list
#' @export

add_partitions <- function(partitions_obj=NULL, to_add=NULL){
	
	ans <- create_partitions_obj(to_add)

	clusterings <- merge(partitions_obj, ans, by=0, all=T, sort=F)
	rownames(clusterings) <- clusterings[, 1]
	clusterings[, 1] <- NULL
	
	return(clusterings)
	
}