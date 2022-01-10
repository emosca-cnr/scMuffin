#' calculate signatures
#' @param genes_by_cells genes_by_cells expression matrix. This matrix must contain positive values only.
#' @param signatures list of gene signatures
#' @param mc.cores number of cores
#' @param nbins number of bins to split the distribution of average gene expression
#' @param nmark_min numner of minimum markers that are required for the succesful calculation of a signature
#' @param ncells_min numner of minimum cells in which a gene set has to be succesfully calculated
#' @param k number of permutations
#' @param kmin minimum number of permutations; due to missing values it is hard to ensure that a signature can be compared to k permutations in every cell
#' @param score_type type of score. if "relative", than the score is the difference between the observed gene set average expression and that of a k permutations; if "mean" the score is equal to the observed gene set average expression
#' @param null_model TRUE if permutations have to be used. Required for score_type="relative"
#' @param verbose verbosity
#' 
#' @export
#' @import Seurat parallel

calculate_signatures <- function(genes_by_cells, signatures=NULL, mc.cores=2, nbins=25, nmark_min = 5, ncells_min = 5, k=100, kmin=50, score_type=c("relative", "mean"), null_model=TRUE, verbose=TRUE){
	
	score_type <- score_type[1]
	#	if(!is.null(custom_signatures)){
	#		signatures <- c(signatures, custom_signatures)
	#	}
	#names(signatures) <- paste0("SIG_", names(signatures))
	print(names(signatures))
	cat("# of signatures: ", length(signatures), "\n")
	
	
	#dataset bins
	data_bins <- NULL
	if(null_model){
		cat("defining bins...\n")
		data_bins <- sc_data_bin(as.matrix(genes_by_cells), nbins = nbins, use.log = TRUE)
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
		
		res_signatures <- lapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin, verbose=verbose))
		
	}else{
		
		res_signatures <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin, verbose=verbose), mc.cores = mc.cores)
		
	}
	
	
	cat("assembling signatures_by_cells matrix...")
	lapply(res_signatures, function(x) print(head(x)))
	
	if(score_type == "relative"){
		SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$avg_delta_score, dimnames = list(rownames(x))))))
	}else{
		SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$case, dimnames = list(rownames(x))))))
	}
	
	return(list(signatures_by_cells=SC_signatures_by_cell_matrix, full=res_signatures, signatures=signatures))
	
}