#' CNV
#' @param x list of dataframes retrieved by 'preprocess_object_for_cnv'.
#' @description Function to be used in calculate_CNV. 
#' @references "Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma" Patel et al., Supplementary materials
#' @author Valentina Nale


CNV <- function(x, wnd_size=100, na.rm=FALSE, genes=NULL) {
	
	N <- length(x)
	
	Ek <- numeric()
	
	wnd_half <- wnd_size/2	
	
	# first check: length
	if (N<(wnd_size+1)) {
		stop('The number of gene is not sufficient to compute CNV')
	}
	
	symbol_names <- genes[(wnd_half+1):(N-(wnd_half))]
	if(na.rm){
		x[x==0] <- NA
	}
	for (i in (wnd_half+1):(N-(wnd_half))) {
		val <- mean(x[(i-wnd_half):(i+wnd_half)], na.rm=na.rm)
		Ek <- c(Ek, val)
		#			Ek <<- Ek
	}
	
	names(Ek)<-symbol_names
	Ek[is.nan(Ek)] <- 0
	
	return(Ek)
}

