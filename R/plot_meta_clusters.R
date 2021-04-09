#' Plot meta clusters
#' @export
#' @import ComplexHeatmap grDevices dendsort

plot_meta_clusters <- function(ov_mat=NULL, meta_clusters=NULL, out_dir="./", file="/meta_clusters.jpg", ...){
	
	
	if(!dir.exists(out_dir)){
		dir.create(out_dir)
	}
	
	grDevices::jpeg(paste0(out_dir, file), width = 180, height = 180, units="mm", res=300)
	
	heat <- ComplexHeatmap::Heatmap(ov_mat, cluster_columns = dendsort::dendsort(meta_clusters$hclust_out), cluster_rows = dendsort::dendsort(meta_clusters$hclust_out), ...)
	draw(heat)
	
	dev.off()
	
	
	
}