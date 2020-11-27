#' calculate_signatures_clusters
#'
#'
#'

calculate_signatures_clusters <- function(sign_list, cell_clusters=cell_clusters, ncells_min = 5, null_model = TRUE){
	
	
	cat("Clusters...\n")
	print(table(cell_clusters))
	
	res_signatures_clusters <- lapply(sign_list, function(i_marker_res) gene_set_score_in_clusters(i_marker_res, cell_clusters=cell_clusters, ncells_min = ncells_min, null_model = null_model))
	
	#signatures-by-clusters matrix
	SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
	
	return(list(full=res_signatures_clusters, signatures_by_clusters=SC_signatures_by_cluster_matrix))
	
}