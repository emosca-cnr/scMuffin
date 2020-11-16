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

boxplot_cluster <- function(features_by_cells, cell_clusters, top_features=NULL, dir_out="./"){
	
	
#  #t-test: for each feature in each cluster
# 	mat <- matrix(0, ncol = length(levels(cell_clusters)), nrow = nrow(features_by_cells), dimnames = list(rownames(features_by_cells), levels(cell_clusters)))
# 	for(i in 1:nrow(mat)){
# 		for(j in 1:ncol(mat)){
# 			mat[i, j] <- t.test(features_by_cells[i, colnames(features_by_cells) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]], features_by_cells[i, !colnames(features_by_cells) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]], alternative = "g")$p.value
# 		}
# 	}
# 	
# 	mat <- apply(mat, 2, p.adjust, method="fdr")
# 	
# 	#top 3 features of each cluster
# 	top_clusters <- lapply(split(data.frame(t(mat)), colnames(mat)), function(x) rownames(mat)[rank(x) <= n_top])
	
	#dir.create(dir_out, showWarnings = F, recursive = T)
	#write.table(mat, file = paste0(dir_out, "/cluster_quantitative_stats.txt"), row.names = T, col,col.names = NA, sep="\t")
	
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
		boxplot(data_clust, at = data_clust_at , las=2, main = paste0("Cluster ", cell_clusters_set[cl]), xaxt = "n", pch=16, outline = F, ylim=c(min(features_by_cells), max(features_by_cells)))
		for(i in 1:nrow(features_by_cells)){
			col <- "pink"
			if(!is.null(top_features)){
				if(rownames(features_by_cells)[i] %in% unlist(top_features[names(top_features) == cell_clusters_set[cl]])){
					col <- "red"
				}
			}
			
			points(jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), data_clust[, i], col=adjustcolor(col, 0.4), pch=16, cex=0.3)
		}
		
		#boxplots of other clusters
		data_no_clust_at <- seq(2, nrow(features_by_cells)*2, by = 2)
		boxplot(data_no_clust, at = seq(2, nrow(features_by_cells)*2, by = 2), las=2, main = paste0("Cluster ", cell_clusters_set[cl]), xaxt = "n", pch=16, add = T, outline = F)
		for(i in 1:nrow(features_by_cells)){
			points(jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), data_no_clust[, i], col=adjustcolor("darkgray", 0.4), pch=16, cex=0.3)
		}
		
		axis(1, data_clust_at, rownames(features_by_cells), las=2, cex=0.5)
	
		dev.off()
	}
}
