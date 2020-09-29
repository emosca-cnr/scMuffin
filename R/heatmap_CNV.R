#' heatmap_CNV
#' @param chr_merged 
#' @param list_ncol_chr
#' @description Plot as an heatmap the resulting data from CNV. 
#' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' @usage heatmap_CNV(chr_merged)
#' @author Valentina Nale
#' @import RColorBrewer Seurat

heatmap_CNV <- function(chr_merged, ngenes_chrom, file="CNV_heatmap.jpg", pal=NULL, n_colors=11, cluster_cells=F, scale_cells=T) {
	
	rotate <- function(x) t(apply(x, 2, rev)) # rotate +90
	
	if(is.null(pal)){
		colors_ <- rev(brewer.pal(n_colors, "RdYlBu"))
	}
		
	ngenes_chrom_cumsum <- cumsum(ngenes_chrom)
	
	if(scale_cells){
		temp <- apply(chr_merged, 2, scale)
		rownames(temp) <- rownames(chr_merged)
		chr_merged <- temp
		rm(temp)
	}
	
	if(cluster_cells){
		seu_obj <- CreateSeuratObject(counts = chr_merged, min.cells = 0, min.features = 0)
		all.genes <- rownames(seu_obj)
		
		seu_obj <- ScaleData(seu_obj, features = all.genes)
		seu_obj <- RunPCA(seu_obj, features = all.genes)
		seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
		seu_obj <- FindClusters(seu_obj)
		
		clusters <- seu_obj@active.ident
		rm(seu_obj)
		
		chr_merged <- chr_merged[, order(clusters)]
	}
	X <- rotate(chr_merged)
	#X <- log2(X)
	
	jpeg(file, width=180, height=180, res=300, units="mm")
	
	layout.show(layout(matrix(c(1, 2), byrow = T, nrow = 1), widths = c(0.8, 0.2))) 
	
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
	abline(h=seq(U[2], U[1], length.out = nrow(chr_merged)+1)[ngenes_chrom_cumsum+1])
	
	
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
