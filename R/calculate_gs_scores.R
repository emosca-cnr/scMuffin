#' Calculate gene set scores
#' @description Calculate gene set scoring
#' @param genes_by_cells genes_by_cells expression matrix. This matrix must contain positive values only.
#' @param gs_list list of gene sets
#' @param mc.cores number of cores
#' @param nbins number of bins to split the distribution of average gene expression
#' @param nmark_min numner of minimum markers that are required for the succesful calculation of a gene set score
#' @param ncells_min numner of minimum cells in which a gene set has to be succesfully calculated
#' @param k number of permutations
#' @param kmin minimum number of permutations; due to missing values it is hard to ensure that a gene set score can be compared to k permutations in every cell
#' @param score_type type of score. if "relative", than the score is the difference between the observed gene set average expression and that of a k permutations; if "mean" the score is equal to the observed gene set average expression
#' @param null_model TRUE if permutations have to be used. Required for score_type="relative"
#' @param verbose verbosity
<<<<<<< HEAD
#' @param na.rm whether to use NA or not
=======
#' @return list of objects containing also the feature object (feat_obj) for visualization purpose
>>>>>>> 0f263d3740a8aecd808cb4e68d2dd82dc5f4fa3b
#' 
#' @export
#' @import Seurat parallel

calculate_gs_scores <- function(genes_by_cells, gs_list=NULL, mc.cores=2, nbins=25, nmark_min = 5, ncells_min = 5, k=100, kmin=50, score_type=c("relative", "mean"), null_model=TRUE, verbose=TRUE, na.rm=TRUE){
	
	score_type <- score_type[1]

	print(names(gs_list))
	cat("# of gene_sets: ", length(gs_list), "\n")
	
	#dataset bins 
	data_bins <- NULL
	if(null_model){
		cat("defining bins...\n")
		data_bins <- sc_data_bin(genes_by_cells, nbins = nbins, na.rm=na.rm)
	}
	
	#signatures-by-cell matrix
	cat("calculating gene set scores...\n")
	if(mc.cores==1){
		
		#debug...
		#res_signatures <- vector("list", length(signatures))
		#for(i in 1:length(signatures)){
		#		cat(i)
		#		res_signatures[[i]] <- gene_set_score(signatures[[i]], genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin)
		#	}
		
	  gs_scores <- lapply(gs_list, function(i_marker_set) gs_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin, verbose=verbose, na.rm=na.rm))
		
	}else{
		
	  gs_scores <- parallel::mclapply(gs_list, function(i_marker_set) gs_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin, verbose=verbose, na.rm=na.rm), mc.cores = mc.cores)
		
	}
	
	
	cat("assembling gene_sets_by_cells matrix...")
	#lapply(gene_set_scores, function(x) print(head(x)))
	
	if(score_type == "relative"){
		gss_by_cell_matrix <- t(do.call(cbind, lapply(gs_scores, function(x) array(x$avg_delta_score, dimnames = list(rownames(x))))))
	}else{
	  gss_by_cell_matrix <- t(do.call(cbind, lapply(gs_scores, function(x) array(x$case, dimnames = list(rownames(x))))))
	}
	
	feat_obj <- create_features_obj(as.data.frame(t(gss_by_cell_matrix)))
	
	return(list(gss_by_cells=gss_by_cell_matrix, by_gs=gs_scores, gs_list=gs_list, feat_obj=feat_obj))
	
}
