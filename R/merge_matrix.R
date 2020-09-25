#' merge_matrix
#' @param features_by_cells features-by-cells input matrix
#' @param expr_score expression scores (output of exp_Rate)
#' @param landscent_list list with the following elements: SR, DPT, potency_states, complete_output (output of landSCENT)
#' @return matrix_Scaled matrix merged and scaled


merge_matrix <- function(features_by_cells, expr_score, landscent_list){
  
  #input
  SR <- landscent_list$SR
  dpt <- landscent_list$DPT
  
  #match colnames
  expr_score <- expr_score[,match(colnames(features_by_cells), colnames(expr_score))]
  SR <- SR[,match(colnames(features_by_cells), colnames(SR))]
  dpt <- dpt[,match(colnames(features_by_cells), colnames(expr_score))]
  
  #merge
  matrix_merged <- rbind(features_by_cells, expR, SR, dpt)
  
  #scale data
  mm_seurat<- CreateSeuratObject(counts = matrix_merged, min.cells = 0, min.features = 0)
  all.genes <- rownames(mm_seurat)
  mm_seurat <- ScaleData(mm_seurat, features = all.genes)
  matrix_Scaled <-GetAssayData(object = mm_seurat, slot = "scale.data")
  
  return(matrix_Scaled)
}