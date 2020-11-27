#' cell cycle - merge the two Tirosh signatures for proliferation analysis
#'
#' @export

proliferation_analysis <- function(genes_by_cells, mc.cores=2, nbins=25, nmark_min = 5, ncells_min = 5, k=100, kmin=50, score_type=c("relative", "mean"), null_model=TRUE){
	
	score_type <- score_type[1]
	
	tirosh <- calculate_signatures(genes_by_cells, signatures=SIG_Tirosh, mc.cores=mc.cores, nbins=nbins, nmark_min = nmark_min, ncells_min = ncells_min, k=k, kmin=kmin, score_type=score_type, null_model=null_model)
	
	ans <- tirosh$signatures_by_cells
	
	if(score_type == "mean"){
		ans <- t(apply(ans, 1, scale))
		colnames(ans) <- colnames(tirosh$signatures_by_cells)
	}
	
	ans <- apply(tirosh$signatures_by_cells, 2, max, na.rm=TRUE) #max between the two signatures
	names(ans) <- colnames(tirosh$signatures_by_cells)
	
	ans[ans==-Inf] <- NA
	
	return(ans)
	
}