#' Calculate CNV
#' @param genes_by_cells genes-by-cells data.frame 
#' @param reference one-column data.frame with a reference expression profile; rownames must match those of genes_by_cell;
#' @param mc.cores number of cores;
#' @param wnd_size number of adjacent genes considered;
#' @param min_genes minimun number of genes expressed in a cell;
#' @param min_cells minimum number of cells in which a gene must be expressed;
#' @param expr_lim min and max values of relative expression; by default, lower and higher whiskers returned by grDevices::boxplot.stats will be used. Set to NA or anything other value such that length(expr_lim) != 2 to disbale the removal of outliers.
#' @param center_cells whether to center cells
#' @param na.rm whether to remove 0 values in CNV estimation
#' @param center_genes whether to center genes or not
#' @param gene_ann optional data.frame with gene annotation with mandatory columns "Chromosome", "symbol" and "pos" (genomic location). If NULL gene locations will be collected from org.Hs.eg.db
#' @importFrom grDevices boxplot.stats
#' @importFrom Matrix Matrix rowSums colSums
#' @description Calculate CNV using the moving average approach firstly described in Patel et al., 2014 Science (DOI: 10.1126/science)
#' @export
#' 
calculate_CNV <- function(genes_by_cells, reference=NULL, mc.cores=1, wnd_size=100, min_genes=200, min_cells=100, expr_lim=NULL, center_cells=TRUE, na.rm=FALSE, center_genes=FALSE, gene_ann=NULL) {
	
 	genes_by_cells[is.na(genes_by_cells)] <- 0
	
 	cat("Genes-by-cells matrix size", dim(genes_by_cells), "\n")
 	
 	cat("Filtering...")
	idx_keep <- rowSums(genes_by_cells !=0) >= min_cells #genes not missing in at least...
	#print(table(idx_keep))
	genes_by_cells <- genes_by_cells[idx_keep, ]

	if(nrow(genes_by_cells) < min_genes){
	  stop("Less than ", min_genes, "available: considering reducing the parameter min_genes.\n")
	}
	
	idx_keep <- colSums(genes_by_cells !=0) >= min(min_genes, nrow(genes_by_cells))
	#print(table(idx_keep))
	genes_by_cells <- genes_by_cells[, idx_keep]
	cat(" done\n")
	cat("Genes-by-cells matrix size", dim(genes_by_cells), "\n")
	
	# ****************************************************************
	if (!is.null(reference)) {
		cat("Adding reference vector...\n")
		
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
		
		if(is.null(expr_lim)){
		  expr_lim <- boxplot.stats(as.numeric(genes_by_cells))$stats[c(1, 5)]
		}
		
		#constraints over relative expression Tirosh et al.
		if(length(expr_lim)==2){
		  cat("applying relative expression limits", expr_lim,"\n")
		  cat("min:", min(genes_by_cells), ", max:", max(genes_by_cells), "\n")
		  genes_by_cells[genes_by_cells < expr_lim[1]] <- expr_lim[1]
		  genes_by_cells[genes_by_cells > expr_lim[2]] <- expr_lim[2]
		  cat("min:", min(genes_by_cells), ", max:", max(genes_by_cells), "\n")
		}
		
	}
	
	
	# ***************************************************************
	
	cat("Preprocess object to calculate CNV...\n")	
	ans <- preprocess_object_for_CNV(genes_by_cells = genes_by_cells, gene_ann=gene_ann)
	
	temp <- factor(names(ans), levels = c(1, 2, 3, 4 ,5 ,6 ,7 ,8 , 9, 10, 11, 12 , 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"))
	temp <- as.character(sort(temp))
	
	ans <- ans[match(temp, names(ans))]
	
	#exclude chrom with < wnd_size genes
	chr_size <- unlist(lapply(ans, nrow))
	cat("Genes in chromosomes:\n")
	print(chr_size)
	if(any(chr_size <= wnd_size)){
	  cat("Excluding chromosomes", names(chr_size)[chr_size<=wnd_size], "\n")
	  cat("Consider reducing wnd_size parameter.\n")
	}
	ans <- ans[chr_size > wnd_size]
	
	#calculate CNV
	cat("Calculating CNV...\n")
	CNV_data_in <- lapply(ans, function(x) Matrix(as.matrix(x), sparse = TRUE))
	
	for(i in 1:length(ans)){
	  
	  cat("Chromosome", names(ans)[i], "...\n")
	  ans[[i]] <- mclapply(setNames(1:ncol(ans[[i]]), colnames(ans[[i]])), function(j_column) CNV(setNames(ans[[i]][, j_column], rownames(ans[[i]])), wnd_size=wnd_size, na.rm = na.rm), mc.cores = mc.cores)
	  ans[[i]] <- do.call(cbind, ans[[i]])
	  
	}
	
	ans <- do.call(rbind, ans)
	
	if(center_cells){
		temp <- apply(ans, 2, scale, scale=F)
		rownames(temp) <- rownames(ans)
		ans <- temp
		rm(temp)
	}
	
	regions_genes <- regions_to_genes(CNV=ans, CNV_input=CNV_data_in)
	
	return(list(CNV=ans, CNV_input=CNV_data_in, regions2genes=regions_genes))
}
