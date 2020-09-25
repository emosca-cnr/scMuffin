#' preprocess_for_heatmap
#' @param result_cnv 
#' @description Bind chromosomes together as a matrix. 
#' @details 
#' @usage preprocess_for_heatmap(result_cnv)
#' @return A matrix to be used by the function heatmap_CNV.
#' @author Valentina Nale

preprocess_for_heatmap <- function(result_cnv) {

	for(i in 1:length(result_cnv)){
		rownames(result_cnv[[i]]) <- paste0("chr", names(result_cnv)[i], "_", rownames(result_cnv[[i]]))
	} 
	
	chr_merged <- do.call(rbind, result_cnv)
	
	return(chr_merged)
}
