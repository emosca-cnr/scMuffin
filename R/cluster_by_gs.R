#' cluster by n-genes
#' @param Seu_obj Seu_obj
#' @param gs vector of genes
#' @import Seurat
#' @export

cluster_by_gs <- function(Seu_obj, gs=NULL){
	
	Seu_obj <- Seurat::ScaleData(Seu_obj, features = gs)
	Seu_obj <- Seurat::RunPCA(Seu_obj, features = gs)
	Seu_obj <- Seurat::FindNeighbors(Seu_obj, dims = 1:10, features = gs)
	Seu_obj <- Seurat::FindClusters(Seu_obj)
	Seu_obj <- Seurat::RunUMAP(Seu_obj, dims = 1:10)
	
	return(Seu_obj)
	
}