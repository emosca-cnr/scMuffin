#' preprocess_for_heatmap
#' @param result_cnv 
#' @description Bind chromosomes together as a matrix. 
#' @details Every chromosome is first transposed. 
#' @usage preprocess_for_heatmap(result_cnv)
#' @return A matrix to be used by the function heatmap_CNV.
#' @author Valentina Nale

heatmap_cnv <- function(result_cnv) {
	list_ncol_chr <- c()
	for (i in 1:length(result_cnv)) {
		temp <- paste0("cnv_chr", i)
		print(temp)
		temp <- result_cnv[[i]]
		temp <- t(temp)	
		ncol_chr <- print(ncol(temp))
		list_ncol_chr <- c(list_ncol_chr, ncol_chr)
		list_ncol_chr <<- list_ncol_chr
		chr_merged <- matrix(nrow = nrow(temp))
		chr_merged <- cbind(chr_merged, temp)
		chr_merged <<- chr_merged
	}
}
