#' Barplot clusters 
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param scMuffinList scMuffinList object
#' @param feature_id names of the features to be used.
#' @param partition_id one among the partitions
#' @param n_features maximum number of features that will be shown
#' @param cex.axis cex.axis
#' @param only_pos_nes whether to consider only positive enrichments
#' @param quantity which column of the hypergeometric result. Default is "p". Other values are "wbd" and "er"
#' @param dir_out string, output directory
#' @importFrom grDevices jpeg
#' @description Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @export

barplot_cluster <- function(scMuffinList=NULL, feature_id=NULL, feature_name=NULL, partition_id=NULL, dir_out="./", n_features=10, only_pos_nes=TRUE, do_scale_features=FALSE, cex.axis=0.8, p.type=c("p", "p_adj")){
	
  if(!dir.exists(dir_out)){
    dir.create(dir_out, recursive = TRUE)
  }
  
  p.type <- match.arg(p.type, c("p", "p_adj"))
  
  cell_category <- scMuffinList[[feature_id]]$summary
  cell_category <- factor(setNames(cell_category[, colnames(cell_category) %in% feature_name], rownames(cell_category)))
  cell_category_names <- setNames(paste0("v", 1:length(levels(cell_category))), levels(cell_category))
  
  cell_clusters <- setNames(scMuffinList$partitions[, partitio_id], rownames(scMuffinList$partitions))
  cell_clusters_set <- levels(cell_clusters)
  
  # clusters-by-feature value table of enrichment p-value
  en_res_p <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "ORA", feature_name = feature_name, quantity = p.type)[[1]]
  en_res_er <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "ORA", feature_name = feature_name, quantity = "er")[[1]]
  en_res_exp <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "ORA", feature_name = feature_name, quantity = "exp")[[1]]
  
  
 	#boxplot for each cluster
	for(cl in 1:length(cell_clusters_set)){ ###for each cluster
		
		grDevices::jpeg(paste0(dir_out, "/cluster_", cell_clusters_set[cl],".jpg"), width=200, height=180, units="mm", res=300)
	  layout(matrix(c(1, 1, 2, 3), nrow = 2), widths = c(0.6, 0.4))
	  par(mar = c(3, 3, 3, 1))
	  par(mgp = c(1.5, .5, 0))
	  
		#distribution of all cells by feature
		#feature data of the cluster
		cell_category_cluster <- table(cell_category,  factor(names(cell_category) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]], levels=c(FALSE, TRUE)))
		cell_category_cluster <- t(cell_category_cluster)
		cell_category_cluster <- rbind(cell_category_cluster, en_res_exp[rownames(en_res_er)==cell_clusters_set[cl], match(colnames(cell_category_cluster), colnames(en_res_er))])
		rownames(cell_category_cluster) <- c("other", "obs", "exp")
		colnames(cell_category_cluster) <- cell_category_names[match(colnames(cell_category_cluster), names(cell_category_names))]
		
		barplot(cell_category_cluster, beside = T, legend.text = T, main=paste("Cluster", cell_clusters_set[cl]), col=pals::alphabet2(3))	

		
		plot(en_res_er[rownames(en_res_er)==cell_clusters_set[cl], ], -log10(en_res_p[rownames(en_res_er)==cell_clusters_set[cl], ]), pch=16, xlab="ER", ylab="-log10(p)", cex.axis=cex.axis)
		plotrix::thigmophobe.labels(en_res_er[rownames(en_res_er)==cell_clusters_set[cl], ], -log10(en_res_p[rownames(en_res_er)==cell_clusters_set[cl], ]), cell_category_names)
		
		plot.new()
		legend(x = "center", legend=paste(cell_category_names, names(cell_category_names)), cex = cex.axis, bty = "n")
		
		
		dev.off()
	}
}
