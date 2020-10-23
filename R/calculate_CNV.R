#' Calculate CNV
#' @param genes_by_cell Preprocessed Seurat object 
#' @param reference_vector NULL. It is possible to add a reference vector, downloaded from GTEx portal.
#' @export
calculate_CNV <- function(genes_by_cells, reference_vector=NULL, mc.cores=2) {
	# ****************************************************************
	if (!is.null(reference_vector)) {
		cat("Adding reference vector...\n")
		genes_by_cells <- GetAssayData(data)
		gbc_mean <- apply(genes_by_cells, 1, mean)
		gtex_mean_final <-apply(reference_vector, 1, mean)
		GTEx_difference <- gtex_mean_final-gbc_mean 
	
		new_a <- as.data.frame(genes_by_cells)
		new_b <- as.data.frame(gtex_mean_final)
	
		merged_data <- merge(new_a,new_b["gtex_mean_final"],by="row.names",all.x=TRUE)
	
		rownames(merged_data) <- merged_data$Row.names
		merged_data$Row.names <- NULL
		genes_by_cells <- merged_data
		
	}
	# ***************************************************************

	cat("Preprocess object to calculate CNV...\n")	
	ans <- preprocess_object_for_CNV(genes_by_cells)
	
	temp <- factor(names(ans), levels = c(1, 2, 3 ,4 ,5 ,6 ,7 ,8 , 9, 10, 11, 12 , 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"))
	temp <- as.character(sort(temp))
	
	ans <- ans[match(temp, names(ans))]
	
	#exclude chrom with <100 genes
	ans <- ans[unlist(lapply(ans, nrow)) > 100]
	
	#calculate CNV
	cat("Calculating CNV...\n")
	ans <- mclapply(ans, function(x) apply(x[, 4:dim(x)[2]], 2, function(y) CNV(y)), mc.cores = mc.cores)

	#merge chromosomes into a unique table
	#ans <- do.call(rbind, ans)
	
	return(ans)
}
