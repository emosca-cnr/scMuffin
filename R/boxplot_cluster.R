#' boxplot_cluster
#'
#' 
#' @importFrom grDevices jpeg
#' @import graphics


boxplot_cluster <- function(features_by_cells_matrix, cell_clusters, n_top=3, dir_out="./"){
	
	
 #t-test: for each features in each cluster
	mat <- matrix(0, ncol = length(levels(cell_clusters)), nrow = nrow(features_by_cells_matrix), dimnames = list(rownames(features_by_cells_matrix), levels(cell_clusters)))
	for(i in 1:nrow(mat)){
		for(j in 1:ncol(mat)){
			mat[i, j] <- t.test(features_by_cells_matrix[i, colnames(features_by_cells_matrix) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]], features_by_cells_matrix[i, !colnames(features_by_cells_matrix) %in% names(cell_clusters)[cell_clusters == colnames(mat)[j]]], alternative = "g")$p.value
		}
	}
	
	mat <- apply(mat, 2, p.adjust, method="fdr")
	
	#top 3 features of each cluster
	top_clusters <- lapply(split(data.frame(t(mat)), colnames(mat)), function(x) rownames(mat)[rank(x) <= n_top])
	
	dir_output <- paste0(dir_out, "/boxplot_cluster/")
	dir.create(dir_output, showWarnings = F)
	
	#boxplot for each cluster
	for(cl in 1:ncol(mat)){
		
		grDevices::jpeg(paste0(dir_output, "cluster_", colnames(mat)[cl],".jpg"), width=180, height=180, units="mm", res=300)
		par(mar = c(10, 4, 2, 1))
		
		#distribution of all cells by feature
		#feature data of the cluster
		data_clust <- as.data.frame(t(features_by_cells_matrix[, colnames(features_by_cells_matrix) %in% names(cell_clusters)[cell_clusters == colnames(mat)[cl]]]))
		#feature data other clusters
		data_no_clust <- as.data.frame(t(features_by_cells_matrix[, ! colnames(features_by_cells_matrix) %in% names(cell_clusters)[cell_clusters == colnames(mat)[cl]]]))

		#boxplots of the cluster
		data_clust_at <- seq(1, nrow(mat)*2, by = 2)
		boxplot(data_clust, at = data_clust_at , las=2, main = paste0("Cluster ", colnames(mat)[cl]), xaxt = "n", pch=16, outline = F, ylim=c(min(features_by_cells_matrix), max(features_by_cells_matrix)))
		for(i in 1:nrow(features_by_cells_matrix)){
			if(rownames(features_by_cells_matrix)[i] %in% top_clusters[[cl]]){
				col <- "red"
			}else{
				col <- "pink"
			}
			points(jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), data_clust[, i], col=adjustcolor(col, 0.4), pch=16, cex=0.3)
		}
		
		#boxplots of other clusters
		data_no_clust_at <- seq(2, nrow(mat)*2, by = 2)
		boxplot(data_no_clust, at = seq(2, nrow(mat)*2, by = 2), las=2, main = paste0("Cluster ", colnames(mat)[cl]), xaxt = "n", pch=16, add = T, outline = F)
		for(i in 1:nrow(features_by_cells_matrix)){
			points(jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), data_no_clust[, i], col=adjustcolor("darkgray", 0.4), pch=16, cex=0.3)
		}
		
		axis(1, data_clust_at, rownames(mat), las=2, cex=0.5)
	
		dev.off()
	}
}
