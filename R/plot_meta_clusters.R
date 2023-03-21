#' Plot meta clusters
#' @param ov_mat overlap matrix
#' @param meta_clusters meta clusters
#' @param out_file output file
#' @param ... further arguments for ComplexHeatmap::Heatmap
#' @export
#' @description Plot meta clusters
#' @import ComplexHeatmap grDevices dendsort

plot_meta_clusters <- function(ov_mat=NULL, meta_clusters=NULL, out_file="meta_clusters.jpg", ...){
	
	
	
	grDevices::jpeg(file, width = 180, height = 180, units="mm", res=300)
	
	heat <- ComplexHeatmap::Heatmap(ov_mat, cluster_columns = dendsort::dendsort(meta_clusters$hclust_out), cluster_rows = dendsort::dendsort(meta_clusters$hclust_out), ...)
	draw(heat)
	
	dev.off()
	
	
	
}