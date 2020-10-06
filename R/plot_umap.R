#' plot_umap
#' @import Seurat

plot_umap <- function(seurat_object, file="umap.jpg", color_by="ident"){

	data_plot <- Seurat::FetchData(genes_by_cells, vars = c("UMAP_1", "UMAP_2", color_by))

	col_levels <- levels(data_plot$data[, 3])
		
	jpeg(file, width=180, height=180, units="mm", res=300)

	par(mar=c(3, 3, 3, 1))
	par(mgp=c(2, 0.7, 0))
	
	layout(matrix(c(1, 2), nrow = 1), widths = c(0.9, 0.1))
	plot(data_plot$data$UMAP_1, data_plot$data$UMAP_2, pch=16, col=rainbow(length(col_levels))[as.numeric(data_plot$data[, 3])], cex=0.4, xlab="UMAP1", ylab = "UMAP2")
	
	par(mar=c(0, 0, 0, 1))
	plot.new()
	legend("center", col_levels, col=rainbow(length(col_levels)), pch=16, bty="n", cex=0.8)
	
	
	dev.off()
	
	
}