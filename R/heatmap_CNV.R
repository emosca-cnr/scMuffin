#' heatmap_CNV
#' @param chr_merged ...
#' @param list_ncol_chr ...
#' @description Plot as an heatmap the resulting data from CNV. 
#' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' @usage heatmap_CNV(chr_merged)
#' @author Valentina Nale
#' @import RColorBrewer Seurat grDevices stats
#' @import graphics
#' @export

heatmap_CNV <- function(chr_merged, ngenes_chrom, file="heatmap_CNV.jpg", pal=NULL, n_colors=11, scale_cells=T, reference=NULL) {
	
	rotate <- function(x) t(apply(x, 2, rev)) # rotate +90
	ans <- NULL
	
	if(is.null(pal)){
		colors_ <- rev(brewer.pal(n_colors, "RdYlBu"))
	}
	
	ngenes_chrom_cumsum <- cumsum(ngenes_chrom)
	
	#scale by column...
	if(scale_cells){
		temp <- apply(chr_merged, 2, scale)
		rownames(temp) <- rownames(chr_merged)
		chr_merged <- temp
		rm(temp)
	}
	
	#clustering
	seu_obj <- CreateSeuratObject(counts = chr_merged, min.cells = 0, min.features = 0)
	all.genes <- rownames(seu_obj)
	
	seu_obj <- ScaleData(seu_obj, features = all.genes)
	seu_obj <- RunPCA(seu_obj, features = all.genes)
	seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
	seu_obj <- FindClusters(seu_obj)
	seu_obj <- BuildClusterTree(seu_obj, features = all.genes, reorder = T)
	hc_cells <- as.hclust(seu_obj@tools$BuildClusterTree)
		
	clusters <- seu_obj@active.ident
	rm(seu_obj)

	ans <- clusters
	
	#order clusters
	clusters_ordered <- clusters[order(clusters)] #re-order columns by cluster
	#order the matrix by clusters
	chr_merged <- chr_merged[, match(names(clusters_ordered), colnames(chr_merged))]
	
	#find the cluster where the reference appears
	if(!is.null(reference) & any(colnames(chr_merged) == reference)){
		ref_cluster <- clusters_ordered[names(clusters_ordered) == reference] #cluster in which the reference occurs
		cat("Reference cluster:", as.character(ref_cluster), "\n")
		cat("Subtracting reference cluster average from CNV profiles...\n")
		
		#update the CNV Matrix, subtracting the average of the reference cluster from CNV profiles
		ref_cluster_avg <- rowMeans(chr_merged[, clusters_ordered==ref_cluster])
		chr_merged <- apply(chr_merged, 2, function(x) x-ref_cluster_avg)
	}
	
	X <- rotate(chr_merged)
	#X <- log2(X)
	
	grDevices::jpeg(file, width=180, height=180, res=300, units="mm")
	
	layout.show(layout(matrix(c(1, 2, 3, 4), byrow = T, nrow = 2), widths = c(0.8, 0.2), heights = c(0.1, 0.9)))
	
	#cluster dendrograms
	par(mar=c(0, 5, 0.1, 1))
	plot(as.dendrogram(hc_cells), type = "rectangle", leaflab = "none", axes=F, edgePar = list())
	
	plot.new()
	
	par(mar=c(5, 5, 1, 1))
	
	if(n_colors==9){
		image(X, xaxt="none", yaxt="none", col=colors_)
	}
	
	if(n_colors==11){
		image(X, xaxt="none", yaxt="none", col=colors_, breaks = c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6)))
	}
	
	
	xx <- seq(0, 1, length.out = nrow(X))
	yy <- seq(0, 1, length.out = ncol(X))
	xx_lab <- rownames(X) #horizontal axis
	yy_lab <- colnames(X) #vertical axis
	#axis(1, xx, xx_lab, las=2, cex=0.3)
	axis(2, rev(yy)[ngenes_chrom_cumsum], gsub("chr([^_]+)_.*", "\\1", rev(yy_lab)[ngenes_chrom_cumsum]), las=2, cex=0.2)
	
	#grid(ncol(X), nrow(X)) #X è la matrice plottata
	U <- par("usr")
	abline(h=seq(U[4], U[3], length.out = nrow(chr_merged)+1)[ngenes_chrom_cumsum+1])
	
	clusters_ordered <- as.numeric(clusters_ordered)
	clust_col_sep <- c(which(clusters_ordered[2:length(clusters_ordered)] - clusters_ordered[1:(length(clusters_ordered)-1)]!=0), length(clusters_ordered)) 
	abline(v=seq(U[1], U[2], length.out = ncol(chr_merged)+1)[clust_col_sep+1]) 	
	
	x_lab <- ans[match(colnames(chr_merged), names(ans))]
	axis(1, xx[clust_col_sep], x_lab[clust_col_sep], cex=0.2)
	
	if(!is.null(reference)){
		ref_cluster_idx <- which(x_lab[clust_col_sep] == ref_cluster) #cluster idx
		
		axis(3, xx[which(rownames(X) == reference)], "R", cex=0.5, font=2)
		axis(1, xx[clust_col_sep][ref_cluster_idx], x_lab[clust_col_sep][ref_cluster_idx], cex=0.2, font=2)
		
	}
	
	# new legend
	par(mar=c(0, 0, 0, 0))
	plot.new()
	
	if(n_colors == 9){
		leg_text <- cbind(seq(min(X), max(X), length.out = n_colors)[-n_colors], seq(min(X), max(X), length.out = n_colors)[-1])
	}
	
	if(n_colors == 11){
		leg_text <- cbind(c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-12], c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-1])
	}
	
	leg_text <- apply(leg_text, 1, function(x) paste0(format(x[1], digits = 2), ", ", format(x[2], digits = 2)))
	legend("center", rev(leg_text), pch=16, col=rev(colors_), cex=0.8, xpd=T, pt.cex=1)
	
	dev.off()
	
	return(clusters)
}
