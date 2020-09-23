#' heatmap_CNV
#' @param chr_merged 
#' @param list_ncol_chr
#' @description Plot as an heatmap the resulting data from CNV. 
#' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' @usage heatmap_CNV(chr_merged)
#' @author Valentina Nale

heatmap_CNV <- function(chr_merged, list_ncol_chr) {
	
	X = chr_merged
	X <- log(X)
	layout.show(layout(matrix(c(1, 2), byrow = T, nrow = 1), widths = c(0.8, 0.2))) 
	
	par(mar=c(5, 5, 1, 1))
	library(RColorBrewer)
	image(X, xaxt="none", yaxt="none", col=rev(brewer.pal(9,"RdYlGn"))) #, breaks = c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6)))
	
	xx <- seq(0, 1, length.out = nrow(X))
	yy <- seq(0, 1, length.out = ncol(X))
	xx_lab <- rownames(X) #horizontal axis
	yy_lab <- colnames(X) #vertical axis
	
	# new legend
	par(mar=c(0, 0, 0, 0))
	plot.new()
	n_colors <- 9
	leg_text <- cbind(seq(min(X), max(X), length.out = n_colors)[-n_colors], seq(min(X), max(X), length.out = n_colors)[-1])
	
	#leg_text <- cbind(c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-12], c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-1])
	leg_text <- apply(leg_text, 1, function(x) paste0(format(x[1], digits = 2), ", ", format(x[2], digits = 2)))
	legend("center", rev(leg_text), pch=16, col=brewer.pal(9, "RdYlGn"), cex=0.8, xpd=T, pt.cex=1)
	
	abline(v=1220, untf = FALSE, col = "blue")
	
	dev.off()
	
	# ****************************************************************************
	# ******************** CHROMOSOME DIVISION / ABLINE (?) **********************
	# ****************************************************************************
	
	
	
	
	
	
	
	
}
