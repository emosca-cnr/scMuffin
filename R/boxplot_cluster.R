#' Bocplt of cluster enrichment in a quantitative feature
#' @param scMuffinList scMuffinList object
#' @param partition_id one among the partitions
#' @param n_features maximum number of features that will be shown
#' @param cex.axis cex.axis
#' @param only_pos_nes whether to consider only positive enrichments
#' @param do_scale_features whether to scale features
#' @param dir_out output directory
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param feature_name the names of the feature that should be considered. It must be one of names(scMuffinList)
#' @description Produce boxplots to visualize the distribution of cell values according to the selected figure. A png figure for each cluster is saved in dir_out.
#' @importFrom grDevices jpeg
#' @importFrom plotrix thigmophobe.labels
#' @import graphics
#' @export

boxplot_cluster <- function(scMuffinList=NULL, feature_name=NULL, partition_id=NULL, dir_out="./", n_features=10, only_pos_nes=TRUE, do_scale_features=FALSE, cex.axis=0.8, width=180, height=180, units="mm", res=300){
	
	if(!dir.exists(dir_out)){
		dir.create(dir_out, recursive = TRUE)
	}
	
  if(length(scMuffinList[[feature_name]]$summary) == 0){
    stop("Can't find scMuffinList[[feature_name]]$summary\n")
  }
  if(!any(colnames(scMuffinList$partitions) == partition_id)){
    stop("Can't find any parition named ", partition_id, "\n")
  }
  
  #cells_by_features <- features$df[, features$type!="factor", drop=F]
  cells_by_features <- scMuffinList[[feature_name]]$summary
  
	cells_by_features <- as.matrix(cells_by_features) ##to fix 
	cells_by_features[is.na(cells_by_features)] <- 0
	cells_by_features <- 	cells_by_features[, colSums(abs(cells_by_features))>0, drop=F]
	
	feature_id <- colnames(cells_by_features)
	
	if(do_scale_features){
		cells_by_features <- scale(cells_by_features)
	}
	
	cell_clusters <- setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))
	cell_clusters_set <- levels(cell_clusters)
	
	
	# clusters-by-feature value table of enrichment p-value
	#en_res_nes <- extract_cluster_enrichment_table(clust_enrich_res, q_type = "nes")
	#en_res_nes <- en_res_nes$cluster_csea_table[names(en_res_nes$cluster_csea_table) == clustering_name][[1]]
	en_res_nes <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "CSEA", feature_id = feature_id, quantity = "nes")
	
	#en_res_q <- extract_cluster_enrichment_table(clust_enrich_res, q_type = "FDRq")
	#en_res_q <- en_res_q$cluster_csea_table[names(en_res_q$cluster_csea_table) == clustering_name][[1]]
	en_res_q <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "CSEA", feature_id = feature_id, quantity = "FDRq")
	
	if(only_pos_nes){
	  en_res_q[en_res_nes < 0] <- 1
	  en_res_nes[en_res_nes < 0] <- 0
	}
	
	#boxplot for each cluster
	ans <- vector("list", length(cell_clusters_set))
	names(ans) <- cell_clusters_set
	
	for(cl in 1:length(cell_clusters_set)){
		#cat(cl)
		
		top_features_ans <- NA
		
		top_features <- list(
			fdr=setNames(en_res_q[rownames(en_res_q)==cell_clusters_set[cl], ], colnames(en_res_q)),
			nes=setNames(en_res_nes[rownames(en_res_nes)==cell_clusters_set[cl], ], colnames(en_res_nes))
		)
		
		#distribution of all cells by feature
		top_features$fdr <- sort(top_features$fdr)
		top_features$nes <- top_features$nes[match(names(top_features$fdr), names(top_features$nes))]

		if(length(top_features$nes)>0){
			
			top_features_ans <- names(top_features$fdr)
			
			top_features <- lapply(top_features, function(x) x[1:min(n_features, length(x))])
			
			cbf_cl <- cells_by_features[, match(names(top_features$fdr), colnames(cells_by_features)), drop=FALSE]
			
			
			grDevices::png(paste0(dir_out, "/cluster_", cell_clusters_set[cl], ".png"), width=width, height=height, units=units, res=res)
			layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = T))
			
			par(mar = c(3, 3, 2, 1))
			par(mgp = c(1.5, .5, 0))
			
			#feature data of the cluster
			cl_cells_idx <- rownames(cbf_cl) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]]
			data_clust <- cbf_cl[cl_cells_idx, , drop=FALSE]
			
			#feature data other clusters
			data_no_clust <- cbf_cl[!cl_cells_idx, , drop=FALSE]
			
			
			data_clust_at <- (seq(2, ncol(cbf_cl)*2, by = 2))
			data_no_clust_at <- (seq(1, ncol(cbf_cl)*2, by = 2))
			
			boxplot_ylim <- c(min(cbf_cl), max(cbf_cl))
			
			#boxplots of the cluster
			box_names <- paste0("f", 1:ncol(data_clust))
			boxplot(data_clust, at = data_clust_at , las=2, main = paste0("Cluster ", cell_clusters_set[cl]), ylab="y", pch=16, outline = F, ylim=boxplot_ylim, xlim=c(0.5, max(data_clust_at)+.5), cex.axis=cex.axis, cex.names=cex.axis, cex.lab=cex.axis, col=NA, names=box_names, cex.main=cex.axis)
			#abline(v=0, lty=2)
			for(i in 1:ncol(cbf_cl)){
				col <- "red"
				points(jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), data_clust[, i], col=adjustcolor(col, 0.8), pch=16, cex=0.3)
			}
			
			#boxplots of other clusters
			boxplot(data_no_clust, at = data_no_clust_at, las=2, main = paste0("Cluster ", cell_clusters_set[cl]), yaxt = "n", pch=16, add = T, outline = F, names=paste0("f", 1:ncol(data_clust), "_ref"), col=NA, cex.axis=cex.axis, cex.names=cex.axis)
			for(i in 1:ncol(cbf_cl)){
				points(jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), data_no_clust[, i] , col=adjustcolor("darkgray", 0.8), pch=16, cex=0.3)
			}
			
			#axis(2, data_clust_at, colnames(data_no_clust), las=2, ...)
			
			#barplot(rev((top_features$nes)), horiz = T, names.arg = "", main = "NES")
			#barplot(rev(-log10(top_features$fdr)), horiz = T, names.arg = "", main = "FDR")
			
			plot(top_features$nes, -log10(top_features$fdr), pch=16, xlab="NES", ylab="-log10(q)", cex.axis=cex.axis, cex.lab=cex.axis, cex=cex.axis)
			plotrix::thigmophobe.labels(top_features$nes, -log10(top_features$fdr), box_names, cex=cex.axis)
			
			plot.new()
			legend(x = "center", legend=paste(box_names, ":", colnames(cbf_cl)), cex = cex.axis, bty = "n", title = "LEGEND")
			
			dev.off()
			
			
		}
		
		ans[[cl]] <- top_features_ans
	}
	
	return(ans)
}
