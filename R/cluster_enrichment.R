#' Calculate cluster enrichment
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features_by_cells matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param dir_out string, output directory
#' @export
#' @author Ettore Mosca

cluster_enrichment <- function(features_by_cells, cell_clusters, n_top=3, dir_out="./"){
	
 #t-test: for each feature in each cluster
	mat <- matrix(0, ncol = length(levels(cell_clusters)), nrow = nrow(features_by_cells), dimnames = list(rownames(features_by_cells), levels(cell_clusters)))
	for(i in 1:nrow(mat)){
		for(j in 1:ncol(mat)){
			mat[i, j] <- t.test(features_by_cells[i, colnames(features_by_cells) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]], features_by_cells[i, !colnames(features_by_cells) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]], alternative = "g")$p.value
		}
	}
	
	mat <- apply(mat, 2, p.adjust, method="fdr")
	
	#GSEA process on features-by-cells
	gsl <- lapply(split(cell_clusters, cell_clusters), function(x) names(x))
	gsea_res <- gsea(t(GetAssayData(features_by_cells, slot="scale.data")), gsl, mc_cores_perm = 2, ord.mode = rep(-1, nrow(features_by_cells)), k = 99)
	
	nes_table <- do.call(cbind, lapply(gsea_res$gs_table, function(x) array(x$nes, dimnames = list(x$id))))
	fdrq_table <- do.call(cbind, lapply(gsea_res$gs_table, function(x) array(x$FDRq, dimnames = list(x$id))))
	
}
