#' heatmap_CNV
#' @param result_cnv 
#' @description Plot as an heatmap the resulting data from CNV. 
#' @details Preliminary part needs to bind the chromosomes, transposed. 
#' @usage heatmap_CNV(result_cnv)
#' @author Valentina Nale

heatmap_cnv <- function(result_cnv) {
	for (i in 1:length(result_cnv)) {
		temp <- paste0("cnv_cnr", i)
		print(temp)
		temp <- result_cnv[[i]]
		temp <- t(temp)	
		chr_merged <- matrix(nrow= nrow(temp))
		chr_merged <- cbind(chr_merged, temp)
		chr_merged <<- chr_merged
	}
}
