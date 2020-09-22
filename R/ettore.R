#' ettore
#' @import parallel
#' @export

ettore <- function(){
	
	# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	# temp <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, keys = head(keys(TxDb.Hsapiens.UCSC.hg19.knownGene)), columns=c('GENEID', 'TXCHROM', 'TXSTART', 'TXEND', 'TXID'), keytype="GENEID")
	# 
	# library(org.Hs.eg.db)
	# gene_locations <- as.data.frame(org.Hs.egCHRLOC)
	# temp <- as.data.frame(org.Hs.egCHRLOCEND)
	# eg2sym <- as.data.frame(org.Hs.egSYMBOL)
	# 
	# gene_locations <- merge(eg2sym, gene_locations, by="gene_id", sort=F)
	# gene_locations <- merge(gene_locations, temp, by=c("gene_id", "Chromosome"), sort=F)
	# gene_locations$pos <- apply(abs(gene_locations[, c("start_location", "end_location")]), 1, min) 
	# gene_locations$start_location <- NULL
	# gene_locations$end_location <- NULL

	
	#expression_clusters <- find_clusters(genes_by_cells)
		
	data_bins <- sc_data_bin(as.matrix(genes_by_cells@assays$RNA@data), nbins = 25, use.log = TRUE)
	
	#	res <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score_in_clusters(i_marker_set, genes_by_cells, expression_clusters, perc=0.25, bins = data_bins, k=100, alt = "two.sided", test="t", nmark_min = 5, ncells_min=10), mc.cores = min(length(signatures), 2))
	
	#res <- mclapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=100, nmark_min = 5, ncells_min = 5), mc.cores = 2)

	res <- lapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=100, nmark_min = 5, ncells_min = 5))
	
	return(res)
	
}