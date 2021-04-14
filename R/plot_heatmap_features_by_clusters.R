#' plot_heatmap_signatures
#' @param signatures_by_clusters list resulting from calculate_signatures_clusters()
#' @param ntop number of top features considered for each cluster
#' @param onlyUp top features are considered only on the basis of their positive deviation from null distribution (up-regulation)
#' @export
#' @import ComplexHeatmap grDevices
#' @importFrom utils write.table


plot_heatmap_features_by_clusters <- function(features_by_clusters, significance_matrix=NULL, sig_threshold=0.05, ntop=10, onlyUp=TRUE, out_dir="./", remove_null_features=FALSE, ...){
	
	
	if(!is.null(significance_matrix)){
		cell_fun_asterisk <- function(j, i, x, y, w, h, fill) {
			if(sig_mat[i, j] < sig_threshold) {
				grid.text("*", x, y)
			}
		}
	}else{
		cell_fun_asterisk <- NULL
	}
	
	
	#output filtering and plotting HEATMAP of each signature type
	if(!is.list(features_by_clusters)){
		features_by_clusters <- list(feature=features_by_clusters)
		if(!is.null(significance_matrix)){
			significance_matrix <- list(feature=significance_matrix)
		}
	}
	for(i in 1:length(features_by_clusters)){
		
		if(is.list(features_by_clusters[[i]])){
			X <- features_by_clusters[[i]]$signatures_by_clusters
		}else{
			X <- features_by_clusters[[i]]
		}
		
		X[is.na(X)] <- 0
		if(remove_null_features){
			X <- X[rowSums(abs(X))>0, ]
			print(dim(X))
		}
		
		
		if(nrow(X)>0){
			
			#find the feature to show in each cluster ###this does not work properly
			top_features <- character()
			for(j in 1:ncol(X)){
				if(onlyUp){
					top_features <- c(top_features, rownames(X)[which(rank(-X[, j]) <= ntop)])
				}else{
					top_features <- c(top_features, rownames(X)[which(rank(-abs(X[, j])) <= ntop)])
				}
			}
			X <- X[rownames(X) %in% top_features, ]
			
			if(!is.null(significance_matrix)){
				sig_mat <- significance_matrix[[i]]
			}
			
			grDevices::jpeg(paste0(out_dir, "/heatmap_", names(features_by_clusters)[i], ".jpg"), width = 180, height = 180, res=300, units="mm")
			
			h_tot_go <- ComplexHeatmap::Heatmap(X, heatmap_legend_param = list(title= names(features_by_clusters)[i]), show_row_names = T, cell_fun = cell_fun_asterisk, ...)
			draw(h_tot_go, heatmap_legend_side = "left")
			
			dev.off()
			
			write.table(X, file=paste0(out_dir, "/table_", names(features_by_clusters)[i], ".txt"), sep = "\t", row.names = T, col.names = NA)
			
			
		}else{
			warning("No available values for", names(features_by_clusters)[i], "\n")
		}
	}
	
}
