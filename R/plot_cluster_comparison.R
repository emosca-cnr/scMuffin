#' plot_cluster_comparison
#' @import Seurat RColorBrewer graphics
#' @export

plot_cluster_comparison <- function(X, dend_row, dend_col, file="heatmap_cluster_comparison.jpg", pal=NULL, n_colors=9, width=180, height=180, res=300, cex.axis=0.6){
	
	rotate <- function(x) t(apply(x, 2, rev)) # rotate +90
	
	if(is.null(pal)){
		colors_ <- brewer.pal(n_colors, "Reds")
	}
	
	X <- rotate(X)
	
	jpeg(file, width = width, height = height, res=res, units="mm")
	
	layout.show(layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = T), heights = c(0.15, 0.95), widths = c(0.2, 0.8)))
	
	par(mar=c(0, 0, 0, 0))
	plot.new()
	par(xpd=T)
	
	leg_text <- cbind(seq(min(X), max(X), length.out = n_colors)[-n_colors], seq(min(X), max(X), length.out = n_colors)[-1])
	leg_text <- apply(leg_text, 1, function(x) paste0(format(x[1], digits = 2), ", ", format(x[2], digits = 2)))
	
	legend("bottom", leg_text, pch=16, col=colors_, cex=0.4, xpd=T, pt.cex=0.6)
	
	c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))
	
	par(mar=c(0, 2, 1, 0))
	plot(dend_col, type = "rectangle", leaflab = "none", axes=F, edgePar = list())
	
	par(mar=c(1.5, 0, 0, 2))
	plot(rev(dend_row), type = "rectangle", leaflab = "none", axes=F, horiz = T)
	#plot.new()
	
	par(mar=c(3, 3, 1.5, 1))
	image(X, xaxt="none", yaxt="none", col=colors_)
	
	xx <- seq(0, 1, length.out = nrow(X))
	yy <- seq(0, 1, length.out = ncol(X))
	xx_lab <- rownames(X) #horizontal axis
	yy_lab <- colnames(X) #vertical axis
	axis(1, xx, xx_lab, las=2, cex.axis=cex.axis)
	axis(2, yy, yy_lab, las=2, cex.axis=cex.axis)
	
	
	for(i in 1:nrow(X)){ #xx
		for(j in 1:ncol(X)){ #yy
			
				text(xx[i], yy[j], format(X[i, j], digits = 2), cex=.6, font = 2)
			
		}
	}
	
	dev.off()
	
	
}