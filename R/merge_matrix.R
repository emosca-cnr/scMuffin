#' merge_matrix
#' 
#' Create and scale a matrix with different features (row) for each cell (column), like signature/landscend scores
#' @param signatures_by_cells matrix, signatures-by-cells input 
#' @param expr_score matrix, expression scores (output of exp_Rate)
#' @param output_landscent dataframe, df with the following elements: SR, DPT, potency_states (output of landSCENT)
#' @return matrix_merged seurat object, object with the merged and scaled matrix
#' @export
#' @author Noemi Di Nanni

merge_matrix <- function(feature_list=NULL){
  
	feature_list <- Reduce(merge, feature_list)
	rownames(feature_list) <- feature_list$id
	feature_list$id <- NULL
	
  return(feature_list)
  
}