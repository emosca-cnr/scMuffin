#' Boxplot clusters 
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features_by_cells matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param dir_out string, output directory
#' @importFrom grDevices jpeg
#' @import graphics
#' @export
#' @author Noemi Di Nanni

boxplot_cluster_ <- function(features, cell_clusters, top_features=NULL, dir_out="./", only_top=FALSE){
	

	features_by_cells <- as.matrix(features$df)
	features_by_cells[is.na(features_by_cells)] <- 0
	features_by_cells <- apply(features_by_cells, 2, scale)
	rownames(features_by_cells) <- rownames(as.matrix(features$df))
	features_by_cells <- t(features_by_cells)
	
	if(only_top){
		features_by_cells <- features_by_cells[rownames(features_by_cells) %in% unique(unlist(lapply(top_features, names))), ]
	}
	
	cell_clusters_set <- levels(cell_clusters)
	
	#boxplot for each cluster
	for(cl in 1:length(cell_clusters_set)){
		
		grDevices::jpeg(paste0(dir_out, "/cluster_", cell_clusters_set[cl],".jpg"), width=180, height=180, units="mm", res=300)
		par(mar = c(10, 4, 2, 1))
		
		#distribution of all cells by feature
		#feature data of the cluster
		data_clust <- as.data.frame(t(features_by_cells[, colnames(features_by_cells) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]]]))
		#feature data other clusters
		data_no_clust <- as.data.frame(t(features_by_cells[, ! colnames(features_by_cells) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]]]))

		#boxplots of the cluster
		data_clust_at <- seq(1, nrow(features_by_cells)*2, by = 2)
		boxplot(data_clust, at = data_clust_at , las=2, main = paste0("Cluster ", cell_clusters_set[cl]), xaxt = "n", pch=16, outline = F, ylim=c(min(features_by_cells), max(features_by_cells)), xlim=c(0.5, max(data_clust_at)+1.5))
		for(i in 1:nrow(features_by_cells)){
			col <- "red"
			if(!is.null(top_features)){
				if(rownames(features_by_cells)[i] %in% unlist(top_features[names(top_features) == cell_clusters_set[cl]])){
					col <- "purple"
				}
			}
			
			points(jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), data_clust[, i], col=adjustcolor(col, 0.8), pch=16, cex=0.3)
		}
		
		#boxplots of other clusters
		data_no_clust_at <- seq(2, nrow(features_by_cells)*2, by = 2)
		boxplot(data_no_clust, at = seq(2, nrow(features_by_cells)*2, by = 2), las=2, main = paste0("Cluster ", cell_clusters_set[cl]), xaxt = "n", pch=16, add = T, outline = F)
		for(i in 1:nrow(features_by_cells)){
			points(jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), data_no_clust[, i], col=adjustcolor("darkgray", 0.8), pch=16, cex=0.3)
		}
		
		axis(1, data_clust_at, rownames(features_by_cells), las=2, cex.axis=0.5)
	
		dev.off()
	}
}
