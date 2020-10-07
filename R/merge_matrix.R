#' merge_matrix
#' @param features_by_cells features-by-cells input matrix
#' @param expr_score expression scores (output of exp_Rate)
#' @param landscent_list list with the following elements: SR, DPT, potency_states, complete_output (output of landSCENT)
#' @return matrix_Scaled matrix merged and scaled


merge_matrix <- function(signatures_by_cells = NULL, expr_score = NULL, cnv = NULL, output_landscent = NULL){
  
	expr_score <- expr_score[, match(colnames(signatures_by_cells), names(expr_score))]
	output_landscent <- t(output_landscent[match(colnames(signatures_by_cells), output_landscent$cell), ])
	
  #merge
  matrix_merged <- rbind(signatures_by_cells, expr_score=expr_score, output_landscent)
  
  #scale data
  matrix_merged <- CreateSeuratObject(counts = matrix_merged, min.cells = 0, min.features = 0)
  all.genes <- rownames(matrix_merged)
  matrix_merged <- ScaleData(matrix_merged, features = all.genes)
  matrix_merged <- GetAssayData(object = matrix_merged, slot = "scale.data")
  
  return(matrix_merged)
  
}