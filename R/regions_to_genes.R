#' Create a map between CNV regions and genes
#' @param CNV CNV matrix
#' @param CNV_input CNV_input
#' @export

regions_to_genes <- function(CNV=NULL, CNV_input=NULL){
	
	###split regions and extract information
	regions <- setNames(strsplit(rownames(CNV), "__"), rownames(CNV))
	
	#get all genes occurring in regions
	check_regions <- unique(unlist(lapply(regions, function(x) x[c(3, 6)])))
	
	### find region-related genes
	genes <- unlist(lapply(CNV_input, rownames))
	genes <- setNames(data.frame(do.call(rbind, strsplit(genes, "__")), stringsAsFactors = F), c("chr", "pos", "symbol"))
	#genes <- unique(unlist(lapply(genes, function(x) x[3])))
	
	#check that all such genes are available in the rows of the initial input
	if(!all(check_regions %in% genes$symbol)){
		warning("Something went wrong with the reconstruction of CNV regions. 'regions2genes' was not calculated. The following genes were not found in the input.\n")
		print(check_regions[!check_regions %in% genes$symbol])
		regions_genes <- NULL
	}else{
		regions_genes <- lapply(regions, function(x) genes[c(which(genes$chr == x[1] & genes$symbol == x[3]) : which(genes$chr == x[4] & genes$symbol == x[6])), ]) #start:end
	}
	
	### assemble the output
	return(regions_genes)
	
}