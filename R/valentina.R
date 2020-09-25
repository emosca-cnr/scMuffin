valentina <- function() {
	punto1 <- preprocess_object_for_CNV(genes_by_cells)
	
	temp <- factor(names(punto1), levels = c(1, 2, 3 ,4 ,5 ,6 ,7 ,8 , 9, 10, 11, 12 , 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"))
	temp <- as.character(sort(temp))
	
	punto1 <- punto1[match(temp, names(punto1))]
	
	#exclude chrom with <100 genes
	punto1 <- punto1[unlist(lapply(punto1, nrow)) > 100]
	
	punto2 <- mclapply(punto1, function(x) apply(x[, 4:dim(x)[2]], 2, function(y) CNV(y)), mc.cores = 2)
	ngenes_chrom <- unlist(lapply(punto2, nrow))
	
	punto3 <- preprocess_for_heatmap(punto2)
	punto4 <- heatmap_CNV(punto3, ngenes_chrom)
	
}
