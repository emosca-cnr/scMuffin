#' plot_umap
#' 
#' Generate a UMAP visualization
#' @param seurat_object seurat object, object with saved dimension reduction components
#' @param file string, file name output
#' @param color_by string, specification of a feature to colour by (e.g. cluster ID)
#' 
#' @import Seurat graphics
#' @export
#' 
plot_umap <- function(seurat_object, file="umap.jpg", labels=NULL, group.by=NULL, feature_plot=FALSE, ...){
	
	jpeg(file, width=180, height=180, units="mm", res=300)
	
	par(mar=c(3, 3, 3, 1))
	par(mgp=c(2, 0.7, 0))
	
	if(feature_plot){
		res <- FeaturePlot(seurat_object, features = group.by, ...)
	}else{
		res <- Seurat::DimPlot(seurat_object, group.by=group.by, ...)
	}
	plot(res)
	
	if(!is.null(labels)){
		
		data_plot <- Seurat::FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", group.by))
		
		labels <- lapply(labels, function(x) paste0(sort(x), collapse = "\n"))
		cluster_xy <- split(data_plot[, 1:2], data_plot[, 3])
		cluster_xy <- do.call(rbind, lapply(cluster_xy, colMeans))
		for(i in 1:nrow(cluster_xy)){
			text(cluster_xy[i, 1], cluster_xy[i, 2], labels[names(labels) == rownames(cluster_xy)[i]][[1]], cex=0.5, font=2)
		}
	}
	
	
	dev.off()
	
	
}