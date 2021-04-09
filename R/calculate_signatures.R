#' calculate signatures
#' 
#' 
#' @export
#' @import Seurat parallel

calculate_signatures <- function(genes_by_cells, signatures=NULL, mc.cores=2, nbins=25, nmark_min = 5, ncells_min = 5, k=100, kmin=50, score_type=c("relative", "mean"), null_model=TRUE){
	
	score_type <- score_type[1]
	#	if(!is.null(custom_signatures)){
	#		signatures <- c(signatures, custom_signatures)
	#	}
	#names(signatures) <- paste0("SIG_", names(signatures))
	print(names(signatures))
	cat("# of signatures: ", length(signatures), "\n")
	
	
	#dataset bins
	data_bins <- sc_data_bin(as.matrix(Seurat::GetAssayData(genes_by_cells)), nbins = nbins, use.log = TRUE)
	
	#signatures-by-cell matrix
	if(mc.cores==1){
		
		#debug...
		#res_signatures <- vector("list", length(signatures))
		#for(i in 1:length(signatures)){
		#		cat(i)
		#		res_signatures[[i]] <- gene_set_score(signatures[[i]], genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin)
		#	}
		
		res_signatures <- lapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin))
		
	}else{
		
		res_signatures <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin), mc.cores = mc.cores)
		
	}
	
	#SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$score_table$avg_delta_score, dimnames = list(rownames(x$score_table))))))
	if(score_type == "relative"){
		SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$avg_delta_score, dimnames = list(rownames(x))))))
	}else{
		SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$case, dimnames = list(rownames(x))))))
	}
	
	return(list(signatures_by_cells=SC_signatures_by_cell_matrix, full=res_signatures))
	
}