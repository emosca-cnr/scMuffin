#' Boxplot clusters 
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param cells_by_features matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param dir_out string, output directory
#' @importFrom grDevices jpeg
#' @import graphics
#' @export
#' @author Noemi Di Nanni

dotplot_cluster <- function(cells_by_features, cell_clusters, top_features=NULL, dir_out="./", cont_tables=NULL){
	
	
	cell_clusters_set <- levels(cell_clusters)
	
	
	#boxplot for each cluster
	for(cl in 1:length(cell_clusters_set)){
		
		grDevices::jpeg(paste0(dir_out, "/cluster_", cell_clusters_set[cl],".jpg"), width=180, height=180, units="mm", res=300)
		par(mar = c(10, 4, 2, 1))
		
		#distribution of all cells by feature
		#feature data of the cluster
		data_clust <- as.data.frame(cells_by_features[rownames(cells_by_features) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]], ])
		
		#feature data other clusters
		data_no_clust <- as.data.frame(cells_by_features[!rownames(cells_by_features) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]], ])
		
		#boxplots of the cluster
		data_clust_at <- seq(1, ncol(cells_by_features)*2, by = 2)
		data_no_clust_at <- seq(2, ncol(cells_by_features)*2, by = 2)
		
		plot(0, xlim=c(1, max(data_no_clust_at)), pch="", ylim = c(1, max(as.numeric(unlist(cells_by_features)))), xaxt="n", xlab="", ylab="", main = paste0("Cluster ", cell_clusters_set[cl]), yaxt="n")
		
		for(i in 1:ncol(cells_by_features)){
			col <- "pink"
			if(!is.null(top_features)){
				if(colnames(cells_by_features)[i] %in% gsub("_.+$", "", unlist(top_features[names(top_features) == cell_clusters_set[cl]]))){
					col <- "red"
				}
			}
			
			data_clust_i <- as.numeric(factor(data_clust[, i], levels = sort(as.numeric(levels(data_clust[, i])))))
			
			points(jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), jitter(data_clust_i, amount = 0.3), col=adjustcolor(col, 0.4), pch=16, cex=0.3)
			j <- (i-1)*length(cell_clusters_set)+cl
			cat(j)
			
			obs_f <- cont_tables[[j]]$observed / sum(cont_tables[[j]]$observed)
			obs_f <- obs_f[order(as.numeric(names(obs_f)))]
			
			text(rep(data_clust_at[i], length(obs_f)), as.numeric(as.factor(names(obs_f))), format(obs_f, digits=2), cex=0.7, font=2)
			text(rep(data_clust_at[i]+0.5, length(obs_f)), as.numeric(as.factor(names(obs_f))), names(obs_f), cex=0.7)
	
		}
		
		
		#dotplots of other clusters
		for(i in 1:ncol(cells_by_features)){
			
			points(jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), jitter(as.numeric(data_no_clust[, i])), col=adjustcolor("darkgray", 0.4), pch=16, cex=0.3)
			j <- (i-1)*length(cell_clusters_set)+cl
			cat(j)
			
			exp_f <- cont_tables[[j]]$expected_orig / sum(cont_tables[[j]]$expected_orig)
			exp_f <- exp_f[order(as.numeric(names(exp_f)))]
			
			text(rep(data_no_clust_at[i], length(exp_f)), as.numeric(as.factor(names(exp_f))), format(exp_f, digits=2), font=2, cex=0.7)
			
		}
		
		axis(1, data_clust_at, colnames(cells_by_features), las=2, cex=0.5)
		
		dev.off()
	}
}
