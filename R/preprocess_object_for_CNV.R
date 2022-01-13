#' preprocess_object_for_CNV
#' @param genes_by_cells genes-by-cells input matrix
#' @description Preprocessing function to obtain a genomically-ordered list of chromosomes. 
#' @details The preliminary step consist of annotation, duplicates and NA values removal. 
#' Then, the matrix is splitted as a list of dataframe, where every dataframe is a chromosome.
#' Chromosomes are ordered from 1 to 22 + X +Y, and then re-ordered by start position. 
#' @usage preprocess_object_for_CNV(genes_by_cells)
#' @return list of genomically-ordered chromosomes
#' @author Valentina Nale
#' @import org.Hs.eg.db
#' @export
preprocess_object_for_CNV <- function(genes_by_cells) {
  
  # retrieve gene informations
  # genes_by_cells=cellObj
  
	cat("Retrieving gene locations...\n")
  gene_locations <- as.data.frame(org.Hs.eg.db::org.Hs.egCHRLOC)
  temp <- as.data.frame(org.Hs.eg.db::org.Hs.egCHRLOCEND)
  eg2sym <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
  
  gene_locations <- merge(eg2sym, gene_locations, by="gene_id", sort=F)
  gene_locations <- merge(gene_locations, temp, by=c("gene_id", "Chromosome"), sort=F)
  gene_locations$pos <- apply(abs(gene_locations[, c("start_location", "end_location")]), 1, min) 
  gene_locations$start_location <- NULL
  gene_locations$end_location <- NULL
  
  matrix_complete <- merge(gene_locations, genes_by_cells, by.x="symbol", by.y=0, sort=FALSE)
  
  # remove chromosome names not in 1:22 and X,Y -- long names, duplicates
  matrix_complete <- matrix_complete[-which(grepl("_", matrix_complete$Chromosome)), ]
  
  # remove NA values across entrezid
  matrix_complete <- matrix_complete[!is.na(matrix_complete$gene_id), ]
  
  # remove duplicates (entregene_id and on row.names('symbol'))
  matrix_reduced <- cbind(matrix_complete, mean = rowMeans(matrix_complete[, 5:dim(matrix_complete)[2]]))
  matrix_reduced <- matrix_reduced[order(-matrix_reduced$mean), ]
  matrix_reduced <- matrix_reduced[!duplicated(matrix_reduced$gene_id),]
  matrix_reduced <- matrix_reduced[!duplicated(matrix_reduced$symbol),]
  matrix_reduced$mean <- NULL
  
  # change rownames with 'SYMBOL'
  mat_rn <- matrix_reduced[,-1]
  rownames(mat_rn) <- matrix_reduced[,1]
  
  # split chromosomes into a list of dataframe
  mat_splitted <- split(mat_rn, mat_rn$Chromosome)
  
  # sort chromosomes by position
  mat_sorted <- lapply(mat_splitted, function(df){
    df[order(df$pos),]
  })
  return(mat_sorted)
}
