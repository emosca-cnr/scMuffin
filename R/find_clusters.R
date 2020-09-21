#' find_clusters
#' @param seurat_object seurat object
#' @param dims number of dimensions
#' @import Seurat
#' @export

find_clusters <- function(seurat_object, dims=1:10){
	
	ans <- Seurat::FindNeighbors(seurat_object, dims = dims)
	ans <- Seurat::FindClusters(ans)
	
	return(ans)
}