#' Calculate signatures clusters
#' @description Calculate signatures clusters
#' @param gs_scores_obj list of gene signature lists
#' @param cell_clusters cell cluster labels
#' @param ncells_min min number of cells required for the clculation of the average signature in the cluster
#' @param null_model TRUE to consider the permutations
#' @return signatures by cluster matrix
#' @export

calculate_gs_scores_in_clusters <- function(gs_scores_obj=NULL, cell_clusters=NULL, ncells_min = 5, null_model = TRUE){
	
	
	cat("Clusters...\n")
	if(!is.factor(cell_clusters)){
		cell_clusters <- as.factor(cell_clusters)
	}
	print(table(cell_clusters))
	
	res_signatures_clusters <- lapply(gs_scores_obj$by_gs, function(i_marker_res) gs_scores_in_clusters(i_marker_res, cell_clusters=cell_clusters, ncells_min = ncells_min, null_model = null_model))
	
	#signatures-by-clusters matrix
	SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
	
	SC_signatures_by_cluster_matrix <- SC_signatures_by_cluster_matrix[, match(levels(cell_clusters), colnames(SC_signatures_by_cluster_matrix))]
	
	return(list(gss_by_clusters=SC_signatures_by_cluster_matrix, by_gs=res_signatures_clusters)) 
	
}