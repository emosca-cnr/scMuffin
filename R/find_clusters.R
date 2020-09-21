#' find_clusters
#' @param seurat_object seurat object
#' @param dims number of dimensions
#' @import Seurat
#' @export

find_clusters <- function(seurat_object){
	
	ans <- Seurat::FindClusters(ans)
	
	return(ans)
}