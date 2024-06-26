#' CNV
#' @param x list of dataframes retrieved by 'preprocess_object_for_cnv'.
#' @param wnd_size number of adjacent genes to consider
#' @param na.rm whether to remove genes that have missing value of expression or not.
#' @description Function to be used in calculate_CNV. 
#' @references "Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma" Patel et al. (DOI: 10.1126/science.1254257), Supplementary materials



CNV <- function(x, wnd_size=100, na.rm=FALSE) {
	
	N <- length(x)
	
	Ek <- numeric()
	Ek_names <- character()
	
	wnd_half <- wnd_size/2	
	
	# first check: length
	if (N<(wnd_size+1)) {
		stop('The number of gene is not sufficient to compute CNV')
	}
	
	#symbol_names <- genes[(wnd_half+1):(N-(wnd_half))]
	if(na.rm){
		x[x==0] <- NA
	}
	for (i in (wnd_half+1):(N-(wnd_half))) {
		val <- mean(x[(i-wnd_half):(i+wnd_half)], na.rm=na.rm)
		Ek <- c(Ek, val)
		Ek_names <- c(Ek_names, paste0(names(x)[c(i-wnd_half, i+wnd_half)], collapse = "__"))
		#			Ek <<- Ek
	}
	
	names(Ek)<-Ek_names
	Ek[is.nan(Ek)] <- 0
	
	#Ek <- matrix(Ek, nrow = 1, dimnames = list(NULL, names(Ek)))
	
	return(Ek)
}

