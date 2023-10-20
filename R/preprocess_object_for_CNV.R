#' preprocess_object_for_CNV
#' @param genes_by_cells genes-by-cells input matrix
#' @param gene_ann optional data.frame with gene annotation with mandatory columns "Chromosome", "symbol" and "pos" (genomic location). If NULL gene locations will be collected from org.Hs.eg.db
#' @description Preprocessing function to obtain a genomically-ordered list of chromosomes. 
#' @details The preliminary step consist of annotation, duplicates and NA values removal. 
#' Then, the matrix is splitted as a list of dataframe, where every dataframe is a chromosome.
#' Chromosomes are ordered from 1 to 22 + X +Y, and then re-ordered by start position. 
#' @return list of genomically-ordered chromosomes
#' @importFrom org.Hs.eg.db org.Hs.egCHRLOC org.Hs.egCHRLOCEND org.Hs.egSYMBOL
#' @export
preprocess_object_for_CNV <- function(genes_by_cells=NULL, gene_ann=NULL) {
  
  # retrieve gene informations
  use_org <- TRUE
  if(!is.null(gene_ann)){
    gene_locations <- gene_ann
    cat("User-provided gene annotations...\n")
    use_org <- FALSE
    if(!all(c("Chromosome", "pos", "symbol") %in% colnames(gene_ann))){
      cat("Can't use the given annotation: proceeding with org.Hs.eg.db...\n")
      use_org <- TRUE
    }
  }
  
  if(use_org){
    cat("Retrieving gene locations...\n")
    gene_locations <- as.data.frame(org.Hs.egCHRLOC)
    temp <- as.data.frame(org.Hs.egCHRLOCEND)
    eg2sym <- as.data.frame(org.Hs.egSYMBOL)
    
    gene_locations <- merge(eg2sym, gene_locations, by="gene_id", sort=F)
    gene_locations <- merge(gene_locations, temp, by=c("gene_id", "Chromosome"), sort=F)
    gene_locations$pos <- apply(abs(gene_locations[, c("start_location", "end_location")]), 1, min) 
    gene_locations$gene_id <- gene_locations$start_location <- gene_locations$end_location <- NULL
    #gene_locations$end_location <- NULL
    gene_locations$Chrom_symbol <- apply(gene_locations[, c("Chromosome", "symbol")], 1, function(x) paste0("chr", x[1], "__", x[2]))
    gene_locations$Chrom_pos_symbol <- apply(gene_locations[, c("Chromosome", "pos", "symbol")], 1, function(x) paste0("chr", x[1], "__", x[2], "__", x[3]))
    gene_locations$Chrom_pos_symbol <- gsub(" +", "", gene_locations$Chrom_pos_symbol )
  }
  
  matrix_complete <- merge(gene_locations, genes_by_cells, by.x="symbol", by.y=0, sort=FALSE)
  
  # remove chromosome names not in 1:22 and X,Y -- long names, duplicates
  matrix_complete <- matrix_complete[-which(grepl("_", matrix_complete$Chromosome)), ]
  
  # remove NA values across symbol
  #matrix_complete <- matrix_complete[!is.na(matrix_complete$gene_id), ]
  matrix_complete <- matrix_complete[!is.na(matrix_complete$symbol), ]
  
  # remove duplicates keep the gene with the highest average expression
  matrix_reduced <- cbind(matrix_complete, mean = rowMeans(matrix_complete[, 6:ncol(matrix_complete)]))
  matrix_reduced <- matrix_reduced[order(-matrix_reduced$mean), ]
  #matrix_reduced <- matrix_reduced[!duplicated(matrix_reduced$gene_id),]
  matrix_reduced <- matrix_reduced[!duplicated(matrix_reduced$Chrom_symbol), ] #it happens that symbol maps to multiple locations
  
  #add rownames and remove columns  
  rownames(matrix_reduced) <- matrix_reduced$Chrom_pos_symbol
  matrix_reduced$symbol <- NULL
  matrix_reduced$Chrom_pos_symbol <- NULL
  matrix_reduced$Chrom_symbol <- NULL
  matrix_reduced$mean <- NULL
  cat("Genes with known location", nrow(matrix_reduced), "\n")
  
  # split chromosomes into a list of dataframe
  matrix_reduced <- split(matrix_reduced, matrix_reduced$Chromosome)
  
  # sort chromosomes by position
  matrix_reduced <- lapply(matrix_reduced, function(x) x[order(x$pos), ])
  matrix_reduced <- lapply(matrix_reduced, function(x) x[, -which(colnames(x) %in% c("Chromosome", "pos"))])
  
  return(matrix_reduced)
}
