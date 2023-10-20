#' Cluster by features
#' 
#' @description Clustering of the features_by_cell matrix
#' @param features_by_cells features by cells
#' @param n_comp numeric, Dimensions of reduction to use as input
#' @param ... arguments passed to Seurat::FindCLusters
#' @return features_by_cells Seurat Object, object with saved dimension reduction components calculate on features by cells matrix
#' @import Seurat
#' @export

cluster_by_features <- function(features_by_cells=NULL, n_comp = 10, ...){
	
	features_by_cells <- CreateSeuratObject(counts = features_by_cells, min.cells = 0, min.features = 0)
	all.genes <- rownames(features_by_cells)
	
	#scale data
	features_by_cells <- ScaleData(features_by_cells, features = all.genes)
	features_by_cells <- RunPCA(features_by_cells, npcs=n_comp, features = all.genes)
	
	n_comp <- min(n_comp, ncol(features_by_cells@reductions$pca))
	
	features_by_cells <- FindNeighbors(features_by_cells, dims = 1:n_comp)
	features_by_cells <- FindClusters(features_by_cells, ...)
	
	# if(cnv){
	# 	#features_by_cells <- Seurat::BuildClusterTree(features_by_cells, features = all.genes, reorder = T)
	# 	#hc_cells <- stats::as.hclust(features_by_cells@tools$BuildClusterTree)
	# 	#ans <- list(clusters=features_by_cells@active.ident, hc=hc_cells, sobj=features_by_cells)
	# 	ans <- list(clusters=features_by_cells@active.ident, sobj=features_by_cells)
	# 	
	# }else{
		ans <- list(clusters=features_by_cells@active.ident, sobj=features_by_cells)
	# }

	return(ans)
	
}
