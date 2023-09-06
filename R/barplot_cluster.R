#' Barplot of cluster enrichment in a categorical feature
#' @param scMuffinList scMuffinList object
#' @param feature_name name of the features to be used.
#' @param feature_id identifier of the feature to be used.
#' @param partition_id one among the partitions
#' @param cex.axis cex.axis
#' @param only_pos_nes whether to consider only positive enrichments
#' @param p.type p for nominal p-value or p_adj for BH FDR.
#' @param dir_out string, output directory
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @importFrom grDevices jpeg
#' @importFrom stats setNames
#' @importFrom plotrix thigmophobe.labels
#' @description Produce barplots (1 for each cluster) of distribution of cells associated with the values of the selected feature. A png figure for each cluster is saved in dir_out.
#' @export

barplot_cluster <- function(scMuffinList=NULL, feature_name=NULL, feature_id=NULL, partition_id=NULL, dir_out="./", only_pos_nes=TRUE, cex.axis=0.8, p.type=c("p", "p_adj"), width=180, height=180, units="mm", res=300){
	
  if(!dir.exists(dir_out)){
    dir.create(dir_out, recursive = TRUE)
  }
  
  p.type <- match.arg(p.type, c("p", "p_adj"))
  
  if(length(scMuffinList[[feature_name]]$summary) == 0){
    stop("Can't find scMuffinList[[feature_name]]$summary\n")
  }
  if(!any(colnames(scMuffinList$partitions) == partition_id)){
    stop("Can't find any parition named ", partition_id, "\n")
  }
  
  cell_category <- scMuffinList[[feature_name]]$summary
  if(!any(colnames(cell_category) %in% feature_id)){
    stop("Can't find any colnames of scMuffinList[[feature_name]]$summary equal to", feature_id, "\n")
  }
  
  cell_category <- factor(setNames(cell_category[, colnames(cell_category) %in% feature_id], rownames(cell_category)))
  cell_category_names <- setNames(paste0("v", 1:length(levels(cell_category))), levels(cell_category))
  
  cell_clusters <- setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))
  cell_clusters_set <- levels(cell_clusters)
  
  # clusters-by-feature value table of enrichment p-value
  en_res_p <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "ORA", feature_id = feature_id, quantity = p.type)[[1]]
  en_res_er <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "ORA", feature_id = feature_id, quantity = "er")[[1]]
  en_res_exp <- extract_cluster_enrichment_table(scMuffinList=scMuffinList, partition_id = partition_id, type = "ORA", feature_id = feature_id, quantity = "exp")[[1]]
  
  
 	#boxplot for each cluster
	for(cl in 1:length(cell_clusters_set)){ ###for each cluster
		
		grDevices::png(paste0(dir_out, "/cluster_", cell_clusters_set[cl],".png"), width=width, height=height, units=units, res=res)
	  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = T))
	  par(mar = c(3, 3, 3, 1))
	  par(mgp = c(1.5, .5, 0))
	  
		#distribution of all cells by feature
		#feature data of the cluster
		cell_category_cluster <- table(cell_category,  factor(names(cell_category) %in% names(cell_clusters)[cell_clusters == cell_clusters_set[cl]], levels=c(FALSE, TRUE)))
		cell_category_cluster <- t(cell_category_cluster)
		cell_category_cluster <- rbind(cell_category_cluster, en_res_exp[rownames(en_res_er)==cell_clusters_set[cl], match(colnames(cell_category_cluster), colnames(en_res_er))])
		rownames(cell_category_cluster) <- c("other", "obs", "exp")
		cell_category_cluster <- cell_category_cluster[-1, ] #remove the other clusters
		cell_category_cluster <- cell_category_cluster[c(2,1), ] #remove the other clusters
		colnames(cell_category_cluster) <- cell_category_names[match(colnames(cell_category_cluster), names(cell_category_names))]
		
		barplot(cell_category_cluster, beside = T, legend.text = T, main=paste("Cluster", cell_clusters_set[cl]), col=pals::alphabet2(2), cex.axis=cex.axis, cex.names=cex.axis, cex.lab=cex.axis, cex.main=cex.axis, args.legend=list(cex=cex.axis))

		
		plot(en_res_er[rownames(en_res_er)==cell_clusters_set[cl], ], -log10(en_res_p[rownames(en_res_er)==cell_clusters_set[cl], ]), pch=16, xlab="ER", ylab="-log10(p)", cex.axis=cex.axis, cex.lab=cex.axis, cex=cex.axis)
		plotrix::thigmophobe.labels(en_res_er[rownames(en_res_er)==cell_clusters_set[cl], ], -log10(en_res_p[rownames(en_res_er)==cell_clusters_set[cl], ]), cell_category_names, cex=cex.axis)
		
		plot.new()
		legend(x = "center", legend=paste(cell_category_names, ":", names(cell_category_names)), cex = cex.axis, bty = "n", title = "LEGEND")
		
		
		dev.off()
	}
}
