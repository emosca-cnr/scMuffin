#' Adds a clustering to the clustering list
#' @param clusterings clusterings object
#' @param clusterings_to_add data.frame of 1 or more cell labels with cell id as row names
#' @export

add_clusterings <- function(clusterings=NULL, clusterings_to_add=NULL){
	
	ans <- create_clusterings(x = clusterings_to_add)

	clusterings <- merge(clusterings, ans, by=0, all=T, sort=F)
	rownames(clusterings) <- clusterings[, 1]
	clusterings[, 1] <- NULL
	
	return(clusterings)
	
}