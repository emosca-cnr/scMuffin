#' preprocess_object_for_CNV
#' @param input_object Seurat object. 
#' @description Preprocessing pipeline to obtain a list of ordered dataframes (one for each chromosome). 
#' @details Preprocessing consist of annotation, removing duplicates and NA values.
#' @usage preprocess_object_for_cnv(input_object)
#' @return list of dataframes, where every dataframe is a chromosome
#' @author Valentina Nale
#' @import Seurat org.Hs.eg.db

preprocess_object_for_CNV <- function(input_object) {
  
  # retrieve gene informations
  # input_object=cellObj
  
  gene_locations <- as.data.frame(org.Hs.egCHRLOC)
  temp <- as.data.frame(org.Hs.egCHRLOCEND)
  eg2sym <- as.data.frame(org.Hs.egSYMBOL)
  
  gene_locations <- merge(eg2sym, gene_locations, by="gene_id", sort=F)
  gene_locations <- merge(gene_locations, temp, by=c("gene_id", "Chromosome"), sort=F)
  gene_locations$pos <- apply(abs(gene_locations[, c("start_location", "end_location")]), 1, min) 
  gene_locations$start_location <- NULL
  gene_locations$end_location <- NULL
  
  # merged Seurat matrix and annotation matrix via "symbol"
  temp_matrix <- Seurat::GetAssayData(object=input_object, slot="data")
  matrix_complete <- merge(gene_locations, temp_matrix, by.x="symbol", by.y=0, sort=FALSE)
  
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
