#' ettore
#' @import parallel
#' @export

ettore <- function(x){
	
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
	
	
	#example data
	load("/home/bioinformatics/ndinanni/SC/SC_2/Yuan2018/PJ016/data_PJ016.RData") #it contains the Seurat object "cell"
	#load("sysdata.rda")
	
	genes_by_cells <- as.matrix(cell@assays$RNA@data)
	
	
	names(seurat_object_ident) <- names(cell$seurat_clusters)
	rm(cell)
	
	res <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score_in_clusters(i_marker_set, 	as.matrix(cell@assays$RNA@data), seurat_object_ident, perc=0.25, bins = data_bins, k=100, alt = "two.sided", test="t", nmark_min = 5, ncells_min=10), mc.cores = min(length(all_markers_sets), 2))
	
	return(res)
	
}