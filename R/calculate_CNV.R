#' Calculate CNV
#' @param genes_by_cells genes-by-cells data.frame 
#' @param reference one-column data.frame with a reference expression profile; rownames must match those of genes_by_cell;
#' @param mc.cores number of cores;
#' @param wnd_size number of adjacent genes considered;
#' @param min_genes minimun number of genes expressed in a cell;
#' @param min_cells minimum numbver of cells in which a gene must be expressed;
#' @param expr_lim min and max values of relative expression; if FALSE this constraints will not be used
#' @param scale_cells whether to scale cells
#' @param na.rm whether to remove 0 values in CNV estimation
#' @param center_genes whether to center genes or not
#' @description Calculate CNV with or without reference vector
#' @export
calculate_CNV <- function(genes_by_cells, reference=NULL, mc.cores=2, wnd_size=100, min_genes=1000, min_cells=100, expr_lim=c(-3, 3), scale_cells=TRUE, na.rm=FALSE, center_genes=FALSE) {
	
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
	#if(center_genes){
	#	gbc_mean <- rowMeans(genes_by_cells)
	#	genes_by_cells <- as.data.frame(apply(genes_by_cells, 2, function(x) x - gbc_mean))
	#}
	
	# ****************************************************************
	if (!is.null(reference)) {
		cat("Adding reference vector...\n")
		
		
		#if(center_genes){
		#	genes_by_cells$gbc_mean <- gbc_mean
		#}
		
		colnames(reference) <- "reference"
		genes_by_cells <- merge(genes_by_cells, reference, by="row.names", sort=F)
		rownames(genes_by_cells) <- genes_by_cells$Row.names
		genes_by_cells$Row.names <- NULL
		
		print("size after merging")
		print(dim(genes_by_cells))
		
		genes_by_cells[is.na(genes_by_cells)] <- 0
		
	}
	
	if(center_genes){
		print("centering genes")
		
		#genes_by_cells$reference <- genes_by_cells$reference - genes_by_cells$gbc_mean
		#genes_by_cells$gbc_mean <- NULL
		
		#gbc_mean <- rowMeans(genes_by_cells)
		#genes_by_cells <- as.data.frame(apply(genes_by_cells, 2, function(x) x - gbc_mean))
		
		temp <- t(apply(genes_by_cells, 1, scale, scale=F))
		colnames(temp) <- colnames(genes_by_cells)
		genes_by_cells <- temp
		rm(temp)
		
	}
	
	#constraints over relative expression Tirosh et al.
	if(length(expr_lim)==2){
		cat("applying relative expression limits\n")
		cat("min:", min(genes_by_cells), ", max:", max(genes_by_cells), "\n")
		genes_by_cells[genes_by_cells < expr_lim[1]] <- expr_lim[1]
		genes_by_cells[genes_by_cells > expr_lim[2]] <- expr_lim[2]
		cat("min:", min(genes_by_cells), ", max:", max(genes_by_cells), "\n")
	}
	
	# ***************************************************************
	
	cat("Preprocess object to calculate CNV...\n")	
	ans <- preprocess_object_for_CNV(genes_by_cells)
	
	temp <- factor(names(ans), levels = c(1, 2, 3 ,4 ,5 ,6 ,7 ,8 , 9, 10, 11, 12 , 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"))
	temp <- as.character(sort(temp))
	
	ans <- ans[match(temp, names(ans))]
	
	#exclude chrom with < wnd_size genes
	ans <- ans[unlist(lapply(ans, nrow)) > wnd_size]
	
	for(i in 1:length(ans)){
		rownames(ans[[i]]) <- paste(rownames(ans[[i]]), ans[[i]]$pos, sep="_")
		ans[[i]] <- ans[[i]][, -c(1:3)]
	}
	
	#calculate CNV
	cat("Calculating CNV...\n")
	ans <- mclapply(ans, function(x) apply(x, 2, function(y) CNV(setNames(y, rownames(x)), wnd_size=wnd_size, na.rm = na.rm)), mc.cores = mc.cores)
	
	#merging chromosomes
	for(i in 1:length(ans)){
		if(!is.matrix(ans[[i]])){
			ans[[i]] <- matrix(ans[[i]], nrow = 1)
		}
		rownames(ans[[i]]) <- paste0("chr", names(ans)[i], "_", rownames(ans[[i]]))
	} 
	
	ans <- do.call(rbind, ans)
	
	if(scale_cells){
		temp <- apply(ans, 2, scale, scale=F)
		rownames(temp) <- rownames(ans)
		ans <- temp
		rm(temp)
	}
	
	return(ans)
}
