#' preprocess_for_heatmap_CNV
#' @param result_cnv object resulting from CNV analyiss
#' @description Bind chromosomes together as a matrix. 
#' @usage preprocess_for_heatmap_CNV(result_cnv)
#' @return A matrix to be used by the function heatmap_CNV.
#' @export

preprocess_for_heatmap_CNV <- function(result_cnv) {

	for(i in 1:length(result_cnv)){
		rownames(result_cnv[[i]]) <- paste0("chr", names(result_cnv)[i], "_", rownames(result_cnv[[i]]))
	} 
	
	chr_merged <- do.call(rbind, result_cnv)
	
	return(chr_merged)
}
