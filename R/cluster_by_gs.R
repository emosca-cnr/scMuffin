#' cluster by n-genes
#' 
#' 

cluster_by_gs <- function(genes_by_cells, gs=NULL){
	
	genes_by_cells <- ScaleData(genes_by_cells, features = gs)
	genes_by_cells <- RunPCA(genes_by_cells, features = gs)
	genes_by_cells <- FindNeighbors(genes_by_cells, dims = 1:10)
	genes_by_cells <- FindClusters(genes_by_cells)
	genes_by_cells <- RunUMAP(genes_by_cells, dims = 1:10)
	
	return(genes_by_cells)
	
}