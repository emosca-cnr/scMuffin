#' merge_matrix
#' 
#' Create and scale a matrix with different features (row) for each cell (column), like signature/landscend scores
#' @param signatures_by_cells matrix, signatures-by-cells input 
#' @param expr_score matrix, expression scores (output of exp_Rate)
#' @param output_landscent dataframe, df with the following elements: SR, DPT, potency_states (output of landSCENT)
#' @return matrix_merged seurat object, object with the merged and scaled matrix
#' @export
#' @author Noemi Di Nanni

merge_matrix <- function(signatures_by_cells = NULL, expr_score = NULL, output_landscent = NULL){
  
	expr_score <- expr_score[match(colnames(signatures_by_cells), names(expr_score))]
	
	if(!is.null(output_landscent)){
		output_landscent <- t(output_landscent[match(colnames(signatures_by_cells), rownames(output_landscent)), c("dpt", "SR")])
	}

  #merge
  matrix_merged <- rbind(signatures_by_cells, expr_score=expr_score, output_landscent)
  
  #scale data
  matrix_merged <- CreateSeuratObject(counts = matrix_merged, min.cells = 0, min.features = 0)
  all.genes <- rownames(matrix_merged)
  matrix_merged <- ScaleData(matrix_merged, features = all.genes)
  #matrix_merged <- GetAssayData(object = matrix_merged, slot = "scale.data")
  
  return(matrix_merged)
  
}