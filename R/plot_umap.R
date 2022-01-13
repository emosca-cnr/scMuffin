#' plot_umap
#' 
#' Generate a UMAP visualization
#' @param Seu_obj seurat object, object with saved dimension reduction components
#' @param labels cluster labels
#' @param group.by a feature to colour by (e.g. cluster ID)
#' @param feature_plot whether to call Seurat::FeaturePlot
#' @param lab_size label size
#' @param lab_color label color
#' @param file string, file name output
#' @param adj_outliers logical, whether to adjust the group.by scores, removing outliers
#' @param ... further arguments for Seurat::FeaturePlot or Seurat::DimPlot
#' @import Seurat graphics ggplot2
#' @export
#' 
plot_umap <- function(Seu_obj, file="umap.jpg", labels=NULL, group.by=NULL, feature_plot=FALSE, lab_size=1, lab_color="black", adj_outliers=FALSE, ...){
	
	if(adj_outliers){
		if(!is.numeric(Seu_obj@meta.data[, colnames(Seu_obj@meta.data) == group.by])){
			message("Cannot adjust ", group.by, "because it's not numeric\n")
		}else{
			Seu_obj@meta.data[, colnames(Seu_obj@meta.data) == group.by] <- adj_outliers_col(Seu_obj@meta.data[, colnames(Seu_obj@meta.data) == group.by])
		}
	}
	
	jpeg(file, width=180, height=180, units="mm", res=300)
	
	par(mar=c(3, 3, 3, 1))
	par(mgp=c(2, 0.7, 0))
	
	if(feature_plot){
		res <- Seurat::FeaturePlot(Seu_obj, features = group.by, ...)
	}else{
		res <- Seurat::DimPlot(Seu_obj, group.by=group.by, ...)
	}
	
	if(!is.null(labels)){
		
		data_plot <- Seurat::FetchData(Seu_obj, vars = c("UMAP_1", "UMAP_2", group.by))
		
		#labels <- lapply(labels, function(x) paste0(sort(x), collapse = "\n")) #remove sort
		labels <- lapply(labels, function(x) paste0(x, collapse = "\n")) #remove sort
		cluster_xy <- split(data_plot[, 1:2], data_plot[, 3])
		cluster_xy <- do.call(rbind, lapply(cluster_xy, colMeans))
		#for(i in 1:nrow(cluster_xy)){
		#	text(cluster_xy[i, 1], cluster_xy[i, 2], labels[names(labels) == rownames(cluster_xy)[i]][[1]], cex=0.5, font=2)
		#}
		plot(res + ggplot2::annotate(geom="text", x=cluster_xy[, 1], y=cluster_xy[, 2], label=labels, color=lab_color, size=lab_size))
	}else{
		plot(res)
	}
	
	
	dev.off()
	
	
}