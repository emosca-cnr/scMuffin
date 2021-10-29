#' Calculate cluster enrichment by chi squared
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features_by_cells matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param fdr fdr threshold, default at 0.05
#' @param dir_out string, output directory
#' @export
#' @author Ettore Mosca

cluster_chisq <- function(cells_by_features, cell_clusters, fdr=0.05, top=2){
	

	mat <- matrix(0, ncol = length(levels(cell_clusters)), nrow = ncol(cells_by_features), dimnames = list(colnames(cells_by_features), levels(cell_clusters)))
	top_states <- mat
	res_tables <- vector("list", nrow(mat)*ncol(mat))
	n<-1
	res_tables_names <- c()
	for(i in 1:nrow(mat)){
		for(j in 1:ncol(mat)){
			
			res_tables_names <- c(res_tables_names, paste0(c(rownames(mat)[i], colnames(mat)[j]), collapse = "_"))
			
			feature_data <- factor(cells_by_features[, i])
			cells_cluster <- rownames(cells_by_features) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]
			cells_other_clusters <- !rownames(cells_by_features) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]
			
			res_tables[[n]] <- table(feature_data[cells_cluster])
			
			table_other <- table(feature_data[cells_other_clusters])
			p_other <- table_other / sum(table_other)
			
			res_tables[[n]] <- chisq.test(res_tables[[n]], p = p_other)
			
			res_tables[[n]]$expected_orig <- table_other
			res_tables[[n]]$contrib <- res_tables[[n]]$residuals^2 / res_tables[[n]]$statistic
			
			mat[i, j] <- res_tables[[n]]$p.value
			top_states[i, j] <- names(res_tables[[n]]$contrib)[which.max(res_tables[[n]]$contrib)]
			n<-n+1
		}
	}
	names(res_tables) <- res_tables_names
	
	if(nrow(mat) > 1){
		mat <- apply(mat, 2, p.adjust, method="fdr")
	}
	
	top_features <- vector("list", ncol(mat))
	for(i in 1:length(top_features)){
		top_features[[i]] <- character()
		temp_top <- which(mat[, i] < fdr & rank(mat[, i] ) <= top)
		if(any(temp_top)){
			for(j in 1:length(temp_top)){
				curr_feat <- rownames(mat)[temp_top[j]]
				top_features[[i]] <- c(top_features[[i]], paste0(c(curr_feat, top_states[temp_top[j], i]), collapse="_"))
			}
		}
		names(top_features) <- colnames(mat)
		
	}
	
	return(list(fdr=mat, ct=res_tables, top_states=top_states, top_features=top_features))
	
}
