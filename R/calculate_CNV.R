#' Calculate CNV
#' @param genes_by_cell Preprocessed Seurat object 
#' @param reference NULL. It is possible to add a reference vector, downloaded from GTEx portal.
#' @export
calculate_CNV <- function(genes_by_cells, reference=NULL, mc.cores=2, wnd_size=100, min_genes=1000, min_cells=100, scale_cells=TRUE, na.rm=FALSE) {
	
	genes_by_cells[is.na(genes_by_cells)] <- 0
	
	cat("Filtering")
	idx_keep <- rowSums(genes_by_cells !=0) >= min_cells #genes not missing in at least...
	print(table(idx_keep))
	genes_by_cells <- genes_by_cells[idx_keep, ]
	
	idx_keep <- colSums(genes_by_cells !=0) >= min(min_genes, nrow(genes_by_cells))
	print(table(idx_keep))
	genes_by_cells <- genes_by_cells[, idx_keep]
	print("size after filtering")
	print(dim(genes_by_cells))
	
	#centering genes-by-cells data
	gbc_mean <- rowMeans(genes_by_cells)
	genes_by_cells <- as.data.frame(apply(genes_by_cells, 2, function(x) x - gbc_mean))
	
	# ****************************************************************
	if (!is.null(reference)) {
		cat("Adding reference vector...\n")
		
		genes_by_cells$gbc_mean <- gbc_mean
		
		colnames(reference) <- "reference"
		genes_by_cells <- merge(genes_by_cells, reference, by="row.names", sort=F)
		print("size after merging")
		print(dim(genes_by_cells))
		
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
	ans <- mclapply(ans, function(x) apply(x[, 4:dim(x)[2]], 2, function(y) CNV(y, wnd_size=wnd_size, genes=rownames(x), na.rm = na.rm)), mc.cores = mc.cores)

	#merging chromosomes
	for(i in 1:length(ans)){
		#ans[[i]] <- as.data.frame(ans[[i]])
		rownames(ans[[i]]) <- paste0("chr", names(ans)[i], "_", rownames(ans[[i]]))
	} 
	
	ans <- do.call(rbind, ans)
	
	if(scale_cells){
		temp <- apply(ans, 2, scale)
		rownames(temp) <- rownames(ans)
		ans <- temp
		rm(temp)
	}
	
	return(ans)
}
