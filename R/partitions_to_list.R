#' Transform partitions to lists
#' @param partitions partitions object
#' @description Clusterings to a list
#' @return List of partitions
#' @export

partitions_to_list <- function(partitions){
	
	#####features-by-clustering
	ans <- split(t(partitions), colnames(partitions))

	for(i in 1:length(ans)){
		ans[[i]] <- array(paste0(names(ans)[i], "_", ans[[i]]), dimnames = list(rownames(partitions)))
	}
	
	return(ans)

}