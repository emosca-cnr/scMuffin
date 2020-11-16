#' Calculate cluster enrichment by chi squared
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features_by_cells matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param dir_out string, output directory
#' @export
#' @author Ettore Mosca

cluster_chisq <- function(cells_by_features, cell_clusters){
	

	mat <- matrix(0, ncol = length(levels(cell_clusters)), nrow = ncol(cells_by_features), dimnames = list(colnames(cells_by_features), levels(cell_clusters)))
	cont_tables <- vector("list", nrow(mat)*ncol(mat))
	n<-1
	for(i in 1:nrow(mat)){
		for(j in 1:ncol(mat)){
			
			cells_cluster <- rownames(cells_by_features) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]
			cont_tables[[n]] <- table(cells_by_features[cells_cluster, i])
			
			cells_other_clusters <- !rownames(cells_by_features) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]
			cont_tables[[n]] <- rbind(clust=cont_tables[[n]], other=table(cells_by_features[cells_other_clusters, i]))
			
			#cont_tables[[n]] <- t(apply(cont_tables[[n]], 1, function(x) x/sum(x)))
			mat[i, j] <- chisq.test(cont_tables[[n]])$p.value
			n<-n+1
		}
	}
	
	mat <- apply(mat, 2, p.adjust, method="fdr")
	
	return(list(fdr=mat, ct=cont_tables))
	
}
