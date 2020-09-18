
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
temp <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, keys = head(keys(TxDb.Hsapiens.UCSC.hg19.knownGene)), columns=c('GENEID', 'TXCHROM', 'TXSTART', 'TXEND', 'TXID'), keytype="GENEID")

library(org.Hs.eg.db)
gene_locations <- as.data.frame(org.Hs.egCHRLOC)
temp <- as.data.frame(org.Hs.egCHRLOCEND)
eg2sym <- as.data.frame(org.Hs.egSYMBOL)

gene_locations <- merge(eg2sym, gene_locations, by="gene_id", sort=F)
gene_locations <- merge(gene_locations, temp, by=c("gene_id", "Chromosome"), sort=F)
gene_locations$pos <- apply(abs(gene_locations[, c("start_location", "end_location")]), 1, min) 
gene_locations$start_location <- NULL
gene_locations$end_location <- NULL
