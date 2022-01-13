#' cell cycle - merge the two Tirosh signatures for proliferation analysis
#' @param genes_by_cells genes_by_cells expression matrix
#' @param mc.cores number of cores
#' @param nbins number of bins to split the distribution of average gene expression
#' @param nmark_min numner of minimum markers that are required for the succesful calculation of a signature
#' @param ncells_min numner of minimum cells in which a gene set has to be succesfully calculated
#' @param k number of permutations
#' @param kmin minimum number of permutations; due to missing values it is hard to ensure that a signature can be compared to k permutations in every cell
#' @param score_type type of score. if "relative", than the score is the difference between the observed gene set average expression and that of a k permutations; if "mean" the score is equal to the observed gene set average expression
#' @param null_model TRUE if permutations have to be used. Required for score_type="relative"
#' @param mean_scale whether to scale the values obtained using score_type="mean"
#' @importFrom utils data
#' @export

proliferation_analysis <- function(genes_by_cells, mc.cores=2, nbins=25, nmark_min = 5, ncells_min = 5, k=100, kmin=50, score_type=c("relative", "mean"), null_model=TRUE, mean_scale=TRUE){
	
	data("SIG_Tirosh", envir=environment())
	score_type <- score_type[1]
	
	SIG_Tirosh <- lapply(SIG_Tirosh, function(x) lapply(x, function(y) unique(y[y %in% rownames(genes_by_cells)])))
	
	tirosh <- calculate_signatures(genes_by_cells, signatures=SIG_Tirosh, mc.cores=mc.cores, nbins=nbins, nmark_min = nmark_min, ncells_min = ncells_min, k=k, kmin=kmin, score_type=score_type, null_model=null_model)
	
	ans <- tirosh$signatures_by_cells
	
	if(score_type == "mean" & mean_scale){
		ans <- t(apply(ans, 1, scale))
		colnames(ans) <- colnames(tirosh$signatures_by_cells)
	}
	
	ans <- apply(ans, 2, max, na.rm=TRUE) #max between the two signatures
	names(ans) <- colnames(tirosh$signatures_by_cells)
	
	ans[ans==-Inf] <- NA
	
	return(ans)
	
}