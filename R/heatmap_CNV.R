#' heatmap_CNV
#' @param cnv CNV matrix
#' @param cnv_clustering cell clustering by CNV
#' @param ref_cluster reference cluster
#' @param file output file
#' @param ... arguments passed to ComplexHeatmap::Heatmap
#' @description Plot an heatmap of the CNV. 
#' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' @author Valentina Nale
#' @import ComplexHeatmap grDevices grid
#' @export

heatmap_CNV <- function(cnv, cnv_clustering, ref_cluster=NULL, file="heatmap_CNV.jpg", ...) {
	
	cnv_clustering$clusters <- cnv_clustering$clusters[order(cnv_clustering$clusters)]
	cnv <- cnv[, match(names(cnv_clustering$clusters), colnames(cnv))]
	# 
	# #find the cluster where the reference appears
	# if(!is.null(reference) & any(colnames(cnv) == reference)){
	# 
	# 	ref_cluster <- cnv_clustering$clusters[names(cnv_clustering$clusters) == reference] #cluster in which the reference occurs
	# 	cat("Reference cluster:", as.character(ref_cluster), "\n")
	# 	cat("Subtracting reference cluster average from CNV profiles...\n")
	# 	
	# 	#update the CNV Matrix, subtracting the average of the reference cluster from CNV profiles
	# 	ref_cluster_avg <- rowMeans(cnv[, cnv_clustering$clusters==ref_cluster])
	# 	cnv <- apply(cnv, 2, function(x) x-ref_cluster_avg)
	# }
	
	
	row_chr <- gsub("(chr[^_]+)_.+", "\\1", rownames(cnv))
	row_chr <- factor(row_chr, levels = unique(row_chr))
	row_chr <- split(row_chr, row_chr)
	ngenes_chr <- unlist(lapply(row_chr, length))
	row_splits <- factor(rep(names(ngenes_chr), ngenes_chr), levels = unique(names(ngenes_chr)))
	col_splits <- cnv_clustering$clusters
	clust_color <- rep("black", length(levels(cnv_clustering$clusters)))
	if(!is.null(ref_cluster)){
		clust_color[levels(cnv_clustering$clusters)==ref_cluster] <- "red"
	}
	column_title_gp <- grid::gpar(col = clust_color)
	
	max_abs <- max(abs(cnv), na.rm = T)
	
	grDevices::jpeg(file, width=180, height=180, res=300, units="mm")
	
	temp<- ComplexHeatmap::Heatmap(cnv, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, row_split = row_splits, row_title_rot=0, row_gap = unit(0, "mm"), column_split = col_splits, column_gap = unit(0, "mm"), border=TRUE, column_title_gp=column_title_gp, name="cnv", ...)
	draw(temp)
	
	dev.off()
	
}
