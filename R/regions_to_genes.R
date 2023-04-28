#' Create a map between CNV regions and genes
#' @param CNV CNV matrix
#' @param CNV_input CNV_input
#' @export

regions_to_genes <- function(CNV=NULL, CNV_input=NULL){
	
	###split regions and extract information
	regions <- setNames(strsplit(rownames(CNV), "__"), rownames(CNV))
	
	#regions <- lapply(regions, function(x) gsub("^chr", "", x))
	
	### find region-related genes
	genes <- unlist(lapply(CNV_input, rownames))
	
	regions_genes <- lapply(regions, function(x) genes[c(which(genes == x[2]):which(genes == x[3]))]) #start:end
	regions_genes <- lapply(regions_genes, function(x) setNames(data.frame(do.call(rbind, strsplit(x, "_")), stringsAsFactors = F), c("symbol", "location")))

	### assemble the output
	return(regions_genes)
	
}