#' clusterings_to_list
#' @param clusterings clusterings object
#' @export

clusterings_to_list <- function(clusterings){
	
	#####features-by-clustering
	ans <- split(t(clusterings), colnames(clusterings))

	for(i in 1:length(ans)){
		ans[[i]] <- array(paste0(names(ans)[i], "_", ans[[i]]), dimnames = list(rownames(clusterings)))
	}
	
	return(ans)

}