#' Calculate CNV
#' 
#' @export
calculate_CNV <- function(genes_by_cells, mc.cores=2) {

	cat("Calculating CNV...\n")	
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
