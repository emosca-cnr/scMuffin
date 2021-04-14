#' Calculate cluster enrichment by gsea approach
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features_by_cells matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param dir_out string, output directory
#' @export
#' @author Ettore Mosca

cluster_gsea <- function(features, cell_clusters, min.cells=100){
	

	X <- as.matrix(features$df)
	X[is.na(X)] <- 0
	
	idx_out <- colSums(X!=0) < min.cells
	names_idx_out <- colnames(X)[idx_out]
	
	X <- X[, !idx_out]
	if(nrow(X)<1){
		stop("not enough cells with values\n")
	}
	
	#GSEA process on features-by-cells
	gsl <- lapply(split(cell_clusters, cell_clusters), function(x) names(x))
	gsea_res <- gsea(X, gsl, mc_cores_perm = 2, ord.mode = rep(-1, ncol(X)), k = 99)
	
	nes_table <- do.call(cbind, lapply(gsea_res$gs_table, function(x) array(x$nes, dimnames = list(x$id))))
	if(any(idx_out)){
		nes_table <- cbind(nes_table, matrix(0, nrow = nrow(nes_table), ncol = length(names_idx_out), dimnames = list(rownames(nes_table), names_idx_out)))
	}
	nes_table <- nes_table[, match(colnames(features$df), colnames(nes_table)), drop=F]
	
	fdrq_table <- do.call(cbind, lapply(gsea_res$gs_table, function(x) array(x$FDRq, dimnames = list(x$id))))
	if(any(idx_out)){
		fdrq_table <- cbind(fdrq_table, matrix(1, nrow = nrow(fdrq_table), ncol = length(names_idx_out), dimnames = list(rownames(fdrq_table), names_idx_out)))
	}
	fdrq_table <- fdrq_table[, match(colnames(features$df), colnames(fdrq_table)), drop=F]
	
	return(list(nes=nes_table, fdrq=fdrq_table))
	
}
