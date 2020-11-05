#' monocle
#' 
#' @export
#' @import monocle

monocle_tree <- function(genes_by_cells, root=NULL){
	
	
	gbc_monocle <- as.cell_data_set(genes_by_cells)
	gbc_monocle <- monocle::reduceDimension(gbc_monocle, norm_method = "none")
	gbc_monocle <- monocle::clusterCells(gbc_monocle, reduction_method = "UMAP")
	gbc_monocle <- monocle::learn_graph(gbc_monocle, use_partition = TRUE)
	gbc_monocle <- monocle::order_cells(gbc_monocle, reduction_method = "UMAP", root_cells = root)
	
	return(gbc_monocle)
	
}