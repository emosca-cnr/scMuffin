#' Boxplot clusters 
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features feature list
#' @param cell_clusters cell clusters
#' @param cluster_enrichment result of assess_cluster_enrichment
#' @param only_top numeric, number of features to be shown with a different color and representing the most significative features according to t-test
#' @param criterion "fdr" to sort features by fdr
#' @param fdr_threshold fdr threshold
#' @param only_pos_nes whether to consider only positive enrichments
#' @param do_scale_features whether to scale features
#' @param dir_out output directory
#' @param ... passed to axis
#' @importFrom grDevices jpeg
#' @import graphics
#' @export

boxplot_cluster <- function(features=NULL, cell_clusters=NULL, cluster_enrichment=NULL, dir_out="./", only_top=10, criterion="fdr", fdr_threshold=NULL, only_pos_nes=TRUE, do_scale_features=FALSE, ...){
	
	if(!dir.exists(dir_out)){
		dir.create(dir_out, recursive = TRUE)
	}
	
	cells_by_features <- as.matrix(features$df)
	cells_by_features[is.na(cells_by_features)] <- 0
	cells_by_features <- 	cells_by_features[, colSums(abs(cells_by_features))>0]
	
	if(do_scale_features){
		cells_by_features <- apply(cells_by_features, 2, scale)
		rownames(cells_by_features) <- rownames(as.matrix(features$df))
	}
	
	cell_clusters_set <- levels(cell_clusters)
	
	if(only_pos_nes){
		cluster_enrichment$fdrq[cluster_enrichment$nes < 0] <- 1
		cluster_enrichment$nes[cluster_enrichment$nes < 0] <- 0
	}
	
	#boxplot for each cluster
	ans <- vector("list", length(cell_clusters_set))
	names(ans) <- cell_clusters_set
	
	for(cl in 1:length(cell_clusters_set)){
		#cat(cl)
		
		top_features_ans <- NA
		
		top_features <- list(
			fdr=cluster_enrichment$fdrq[rownames(cluster_enrichment$fdrq)==cell_clusters_set[cl], ],
			nes=cluster_enrichment$nes[rownames(cluster_enrichment$nes)==cell_clusters_set[cl], ]
		)
		
		#distribution of all cells by feature
		if(criterion=="fdr"){
			top_features$fdr <- sort(top_features$fdr)
			top_features$nes <- top_features$nes[match(names(top_features$fdr), names(top_features$nes))]
			if(!is.null(fdr_threshold)){
				idx <- top_features$fdr < fdr_threshold
				top_features$fdr <- top_features$fdr[idx]
				top_features$nes <- top_features$nes[idx]
			}
		}else{
			top_features$nes <- sort(abs(top_features$nes[rownames(top_features$nes)==cell_clusters_set[cl], ]), decreasing = T)
			top_features$fdr <- top_features$fdr[match(names(top_features$nes), names(top_features$fdr))]
		}
		
		if(length(top_features$nes)>0){
			
			top_features_ans <- names(top_features$fdr)
			
			top_features <- lapply(top_features, function(x) x[1:min(only_top, length(x))])
			
			cbf_cl <- cells_by_features[, match(names(top_features$fdr), colnames(cells_by_features)), drop=FALSE]
			
			
			grDevices::jpeg(paste0(dir_out, "/cluster_", cell_clusters_set[cl],".jpg"), width=200, height=180, units="mm", res=300)
			layout(matrix(c(1, 2, 3), nrow = 1, byrow = F), widths = c(0.6, 0.2, 0.2))
			
			par(mar = c(4, 10, 2, 1))
			
			#feature data of the cluster
			cl_cells_idx <- rownames(cbf_cl) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]]
			data_clust <- cbf_cl[cl_cells_idx, , drop=FALSE]
			
			#feature data other clusters
			data_no_clust <- cbf_cl[!cl_cells_idx, , drop=FALSE]
			
			
			data_clust_at <- rev(seq(2, ncol(cbf_cl)*2, by = 2))
			data_no_clust_at <- rev(seq(1, ncol(cbf_cl)*2, by = 2))
			
			#boxplots of the cluster
			boxplot(data_clust, at = data_clust_at , las=1, main = paste0("Cluster ", cell_clusters_set[cl]), yaxt = "n", pch=16, outline = F, ylim=c(min(cbf_cl), max(cbf_cl)), xlim=c(0.5, max(data_clust_at)+.5), horizontal=TRUE)
			abline(v=0, lty=2)
			for(i in 1:ncol(cbf_cl)){
				col <- "red"
				points(data_clust[, i], jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), col=adjustcolor(col, 0.8), pch=16, cex=0.3)
			}
			
			#boxplots of other clusters
			boxplot(data_no_clust, at = data_no_clust_at, las=1, main = paste0("Cluster ", cell_clusters_set[cl]), yaxt = "n", pch=16, add = T, outline = F, horizontal=TRUE)
			for(i in 1:ncol(cbf_cl)){
				points(data_no_clust[, i], jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), col=adjustcolor("darkgray", 0.8), pch=16, cex=0.3)
			}
			
			axis(2, data_clust_at, colnames(data_no_clust), las=2, ...)
			
			par(mar = c(4, 1, 2, 1))
			barplot(rev((top_features$nes)), horiz = T, names.arg = "", main = "NES")
			
			barplot(rev(-log10(top_features$fdr)), horiz = T, names.arg = "", main = "FDR")
			
			dev.off()
			
			
		}
		
		ans[[cl]] <- top_features_ans
	}
	
	return(ans)
}
