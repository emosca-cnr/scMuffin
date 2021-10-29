#' Get cells in each meta cluster
#' @param  cl_list cluterings list
#' @param meta_clusters meta clusters
#' @export
#' @importFrom stats setNames
#' 

get_meta_clusters <- function(cl_list, meta_clusters){

	meta_clusters$clusters$type <- gsub("_[^_]+$", "", meta_clusters$clusters$cluster)
	
	ans <- vector("list", nrow(meta_clusters$clusters))
	
	for(i in 1:nrow(meta_clusters$clusters)){
		
		idx_cl_list <- meta_clusters$clusters$type[i]
		ans[[i]] <- data.frame(meta_cl=meta_clusters$clusters$meta_cl[i], type=meta_clusters$clusters$type[i], cluster_id=meta_clusters$clusters$cluster[i], cell_id=names(cl_list[[idx_cl_list]])[cl_list[[idx_cl_list]] == meta_clusters$clusters$cluster[i]], stringsAsFactors = F)
		
	}
	
	ans <- Reduce(rbind, ans)
	
	temp <- split(ans, ans$meta_cl)
	temp <- lapply(temp, function(x) table(x$cell_id, x$cluster_id))
	ans <- list(clusters=unique(ans[, c("meta_cl", "cell_id")]), occurrence=temp)
	
	return(ans)	
}