#' re_clustering
#' 
#' Clustering of the features_by_cell matrix
#' @param features_by_cells matrix, features by cells matrix 
#' @param n_comp numeric, Dimensions of reduction to use as input
#' @return features_by_cells Seurat Object, object with saved dimension reduction components calculate on features by cells matrix
#' @import Seurat
#' @export
#' @author Noemi Di Nanni

re_clustering <- function(features_by_cells, n_comp = 10){
	
	features_by_cells <- RunPCA(features_by_cells, features = rownames(features_by_cells))

	features_by_cells <- FindNeighbors(features_by_cells, dims = 1:n_comp)
	features_by_cells <- FindClusters(features_by_cells)
	features_by_cells <- RunUMAP(features_by_cells, dims = 1:n_comp)
	
	return(features_by_cells)
	
}
