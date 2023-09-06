#' Define a proliferation score
#' @param scMuffinList scMuffinList object
#' @param mc.cores number of cores
#' @param nbins number of bins to split the distribution of average gene expression
#' @param nmark_min numner of minimum markers that are required for the succesful calculation of a signature
#' @param ncells_min numner of minimum cells in which a gene set has to be succesfully calculated
#' @param k number of permutations
#' @param kmin minimum number of permutations; due to missing values it is hard to ensure that a signature can be compared to k permutations in every cell
#' @param score_type type of score. if "relative", than the score is the difference between the observed gene set average expression and that of a k permutations; if "mean" the score is equal to the observed gene set average expression
#' @param null_model TRUE if permutations have to be used. Required for score_type="relative"
#' @param mean_scale whether to scale the values obtained using score_type="mean"
#' @param gsl gene set list with cell cycle genes
#' @importFrom utils data
#' @description Define a proliferation score on the basis of two cell cycle gene sets. See SIG_Tirosh
#' @return scMuffinList with element "proliferation", a list with:
#' \itemize{
#' \item{summary: data.frame with proliferation value}
#' \item{full: data.frame with the vales of each of the two considered gene sets}
#' }
#' @export

proliferation_analysis <- function(scMuffinList = NULL, mc.cores=1, nbins=25, nmark_min = 5, ncells_min = 5, k=99, kmin=49, score_type=c("relative", "mean"), null_model=TRUE, mean_scale=TRUE, gsl=NULL){
	
	score_type <- score_type[1]
	gsls_EntrezID <- gsls_Symbol <- NULL #to please the check
	
	if(is.null(gsl)){
		cat("Using gsls_Symbol$Tirosh as gene sets for G1S and G2M.\n")
		gsl <- gsls_Symbol$Tirosh
	}
	
	names_gsl <- names(gsl)
	
	gsl <- lapply(gsl, function(x) lapply(x, function(y) unique(y[y %in% rownames(scMuffinList$normalized)])))
	print(lengths(gsl))
	
	scMuffinList <- calculate_gs_scores(scMuffinList, gs_list=gsl, mc.cores=mc.cores, nbins=nbins, nmark_min = nmark_min, ncells_min = ncells_min, k=k, kmin=kmin, score_type=score_type, null_model=null_model)
	
	#ans <- tirosh$gss_by_cells
	#NA -> 0
	scMuffinList$gene_set_scoring$summary[names_gsl[1]][is.na(scMuffinList$gene_set_scoring$summary[names_gsl[1]])] <- 0
	scMuffinList$gene_set_scoring$summary[names_gsl[2]][is.na(scMuffinList$gene_set_scoring$summary[names_gsl[2]])] <- 0
	
	if(score_type == "mean" & mean_scale){
	  scMuffinList$gene_set_scoring$summary[, names_gsl] <- t(scale(t(scMuffinList$gene_set_scoring$summary[, names_gsl])))
	}

	scMuffinList$proliferation <- list(
	  summary=data.frame(Proliferation_score=apply(scMuffinList$gene_set_scoring$summary[, names_gsl], 1, max, na.rm=TRUE)), #max between the two signatures
	  full=scMuffinList$gene_set_scoring$summary[, names_gsl]
	)
	
	scMuffinList$gene_set_scoring$summary[names_gsl[1]] <- scMuffinList$gene_set_scoring$summary[names_gsl[2]] <- NULL
	scMuffinList$gene_set_scoring$full[names_gsl] <- NULL

	#ans[ans==-Inf] <- NA
	
	return(scMuffinList)
	
}
