find_baseline <- function(cnv_score, n_signal=10, n_r=10, signal=NULL, r=NULL){
	
	t_signal <- cnv_score$CNV_signal[order(cnv_score$CNV_signal)[n_signal]]
	t_r <- cnv_score$CNV_R[order(cnv_score$CNV_R)[n_signal]]
	
	ans <- list(
		cells=rownames(cnv_score)[cnv_score$CNV_signal <= t_signal | cnv_score$CNV_R <= t_r],
		t_signal=t_signal,
		t_r=t_r
	)
	
	
	return(ans)
	
}