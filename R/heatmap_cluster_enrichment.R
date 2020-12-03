#' heatmap of signatures
#' @param X ...
#' @param file ...
#' @description Plot as an heatmap the resulting data from CNV. 
#' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' @usage heatmap_CNV(chr_merged)
#' @author Ettore Mosca
#' @import RColorBrewer graphics
## @importFrom dendextend color_branches
#' @export

heatmap_cluster_enrichment <- function(X, Y, file="heatmap_cluster_enrichment.jpg", pal=NULL, n_colors=11, seurat_dendrogram=NULL, width=180, height=180, res=300, cex.axis=0.3, fdr_cutoff=0.05) {
	
	rotate <- function(x) t(apply(x, 2, rev)) # rotate +90
	
	if(is.null(pal)){
		colors_ <- rev(brewer.pal(n_colors, "RdYlBu"))
	}
	

	#cells clustering
	if(!is.null(seurat_dendrogram)){
		hc_col <- seurat_dendrogram
		X <- X[, match(hc_col$labels, gsub("^.+_", "", colnames(X)))]
	}else{
		hc_col <- hclust(dist(t(X))) #COLUMNS
		X <- X[, hc_col$order]
	}
	
	#signatures clustering
	hc_row <- hclust(dist(X)) #ROWS
	X <- X[hc_row$order, ]
	
		
	sample_labels <- apply(X, 2, which.max)
	
	clust_col <- cutree(hc_col, 7)
	clust_col <- clust_col[match(colnames(X), names(clust_col))]
	clust_col_sep <- which(clust_col[2:length(clust_col)] - clust_col[1:(length(clust_col)-1)]!=0)
	#dend_col <- dendextend::color_branches(as.dendrogram(hc_col), k = 7)
	dend_col <- as.dendrogram(hc_col)
	
	clust_row <- cutree(hc_row, 5)
	clust_row <- clust_row[match(rownames(X), names(clust_row))]
	clust_row_sep <- which(clust_row[2:length(clust_row)] - clust_row[1:(length(clust_row)-1)]!=0)
	#dend_row <- dendextend::color_branches(as.dendrogram(hc_row), k = 5)
	dend_col <- as.dendrogram(hc_col)
	
	#rotation
	X <- rotate(X)
	Y <- rotate(Y)
	Y <- Y[match(rownames(X), rownames(Y)), match(colnames(X), colnames(Y))]
	
	jpeg(file, width = width, height = height, res=res, units="mm")
	
	layout.show(layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = T), heights = c(0.15, 0.95), widths = c(0.2, 0.8)))
	
	par(mar=c(0, 0, 0, 0))
	plot.new()
	par(xpd=T)
	
	leg_text <- cbind(c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-12], c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-1])
	leg_text <- apply(leg_text, 1, function(x) paste0(format(x[1], digits = 2), ", ", format(x[2], digits = 2)))
	
	legend("bottom", leg_text, pch=16, col=colors_, cex=0.4, xpd=T, pt.cex=0.6)
	
	c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))
	
	par(mar=c(0, 2, 1, 0))
	plot(dend_col, type = "rectangle", leaflab = "none", axes=F, edgePar = list())
	
	par(mar=c(1.5, 0, 0, 2))
	plot(rev(dend_row), type = "rectangle", leaflab = "none", axes=F, horiz = T)
	#plot.new()
	
	par(mar=c(6, 3, 1.5, 1))
	image(X, xaxt="none", yaxt="none", col=colors_, breaks = c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6)))
	
	xx <- seq(0, 1, length.out = nrow(X))
	yy <- seq(0, 1, length.out = ncol(X))
	xx_lab <- rownames(X) #horizontal axis
	yy_lab <- colnames(X) #vertical axis
	axis(1, xx, xx_lab, las=2, cex.axis=cex.axis)
	axis(2, yy, yy_lab, las=2, cex.axis=cex.axis)
	
	for(i in 1:nrow(X)){ #xx
		for(j in 1:ncol(X)){ #yy
			
			if(Y[i, j] < fdr_cutoff & X[i, j] > 0){ #significant and positive nes
				points(xx[i], yy[j], pch="*", cex=2, lwd=2, col="pink")
			}
			
		}
	}
	
	dev.off()
	
}
