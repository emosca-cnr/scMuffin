#' calculate_signatures_clusters
#' @param sign_list list of gene signature lists
#' @param cell_clusters cell cluster labels
#' @param ncells_min min number of cells required for the clculation of the average signature in the cluster
#' @param null_model TRUE to consider the permutations
#'
#' @export

calculate_signatures_clusters <- function(sign_list, cell_clusters=cell_clusters, ncells_min = 5, null_model = TRUE){
	
	
	cat("Clusters...\n")
	if(!is.factor(cell_clusters)){
		cell_clusters <- as.factor(cell_clusters)
	}
	print(table(cell_clusters))
	
	res_signatures_clusters <- lapply(sign_list, function(i_marker_res) gene_set_score_in_clusters(i_marker_res, cell_clusters=cell_clusters, ncells_min = ncells_min, null_model = null_model))
	
	#signatures-by-clusters matrix
	SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
	
	SC_signatures_by_cluster_matrix <- SC_signatures_by_cluster_matrix[, match(levels(cell_clusters), colnames(SC_signatures_by_cluster_matrix))]
	
	return(list(full=res_signatures_clusters, signatures_by_clusters=SC_signatures_by_cluster_matrix))
	
}