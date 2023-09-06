#' Create a map between CNV regions and genes
#' @param CNV CNV matrix
#' @param CNV_input CNV_input
#' @export

regions_to_genes <- function(CNV=NULL, CNV_input=NULL){
	
	###split regions and extract information
	regions <- setNames(strsplit(rownames(CNV), "__"), rownames(CNV))
	
	#get all genes occurring in regions
	check_regions <- unique(unlist(lapply(regions, function(x) x[-1])))
	
	### find region-related genes
	genes <- unlist(lapply(CNV_input, rownames))
	
	#check that all such genes are available in the rows of the initial input
	if(!all(check_regions %in% genes)){
		warning("Something went wrong with the reconstruction of CNV regions. 'regions2genes' was not calculated. The following genes were not found in the input.\n")
		print(check_regions[!check_regions %in% genes])
		regions_genes <- NULL
	}else{
		regions_genes <- lapply(regions, function(x) genes[c(which(genes == x[2]):which(genes == x[3]))]) #start:end
		regions_genes <- lapply(regions_genes, function(x) setNames(data.frame(do.call(rbind, strsplit(x, "_")), stringsAsFactors = F), c("symbol", "location")))
	}
	
	### assemble the output
	return(regions_genes)
	
}