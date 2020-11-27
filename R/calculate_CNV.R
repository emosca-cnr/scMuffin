#' Calculate CNV
#' @param genes_by_cell Preprocessed Seurat object 
#' @param reference NULL. It is possible to add a reference vector, downloaded from GTEx portal.
#' @export
calculate_CNV <- function(genes_by_cells, reference=NULL, mc.cores=2, wnd_size=100) {
	
	genes_by_cells[is.na(genes_by_cells)] <- 0
	
	#centering genes-by-cells data
	gbc_mean <- rowMeans(genes_by_cells)
	genes_by_cells <- as.data.frame(apply(genes_by_cells, 2, function(x) x - gbc_mean))
	
	# ****************************************************************
	if (!is.null(reference)) {
		cat("Adding reference vector...\n")
		
		genes_by_cells$gbc_mean <- gbc_mean
		
		colnames(reference) <- "reference"
		genes_by_cells <- merge(genes_by_cells, reference, by="row.names", all.x=TRUE, sort=F)
		
		genes_by_cells$reference <- genes_by_cells$reference - genes_by_cells$gbc_mean
		genes_by_cells$gbc_mean <- NULL
		
		rownames(genes_by_cells) <- genes_by_cells$Row.names
		genes_by_cells$Row.names <- NULL
		
		genes_by_cells[is.na(genes_by_cells)] <- 0

	}
	# ***************************************************************

	cat("Preprocess object to calculate CNV...\n")	
	ans <- preprocess_object_for_CNV(genes_by_cells)
	
	temp <- factor(names(ans), levels = c(1, 2, 3 ,4 ,5 ,6 ,7 ,8 , 9, 10, 11, 12 , 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"))
	temp <- as.character(sort(temp))
	
	ans <- ans[match(temp, names(ans))]
	
	#exclude chrom with <100 genes
	ans <- ans[unlist(lapply(ans, nrow)) > wnd_size]
	
	#calculate CNV
	cat("Calculating CNV...\n")
	ans <- mclapply(ans, function(x) apply(x[, 4:dim(x)[2]], 2, function(y) CNV(y, wnd_size=wnd_size)), mc.cores = mc.cores)

	#merge chromosomes into a unique table
	#ans <- do.call(rbind, ans)
	
	return(ans)
}
