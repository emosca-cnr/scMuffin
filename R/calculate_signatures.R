#' calculate signatures
#' 
#' 
#' @export
#' @import Seurat parallel

calculate_signatures <- function(genes_by_cells, custom_signatures=NULL, mc.cores=2, cell_clusters=NULL){
	
	if(!is.null(custom_signatures)){
		signatures <- c(signatures, custom_signatures)
	}
	names(signatures) <- paste0("SIG_", names(signatures))
	cat("# of signatures: ", length(signatures), "\n")
	
	if(is.null(cell_clusters)){
		cell_clusters <- genes_by_cells@active.ident
	}
	cell_clusters <- cell_clusters[match(colnames(genes_by_cells), names(cell_clusters))]
	
	cat("Clusters...\n")
	print(table(cell_clusters))
	
	
	#dataset bins
	data_bins <- sc_data_bin(as.matrix(Seurat::GetAssayData(genes_by_cells)), nbins = 25, use.log = TRUE)
	
	#signatures-by-cell matrix
	res_signatures <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=100, nmark_min = 5, ncells_min = 5), mc.cores = mc.cores)
	
	SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$score_table$avg_delta_score, dimnames = list(rownames(x$score_table))))))
	
	#gene-set score per cluster list
	res_signatures_clusters <- lapply(res_signatures, function(i_marker_res) gene_set_score_in_clusters(i_marker_res$score_table, cell_clusters=cell_clusters, ncells_min = 5))
	
	#signatures-by-clusters matrix
	SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
	
	return(list(signatures_by_cells=SC_signatures_by_cell_matrix, signatures_by_clusters=SC_signatures_by_cluster_matrix))
	
}