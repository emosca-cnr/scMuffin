#' gene_set_score_in_clusters
#' @param matrix matrix of features
#' @param n_comp Dimensions of reduction to use as input
#' @return features_by_cells Seurat Object
#' 
#' @import Seurat
#' @export
#' @Author Noemi Di Nanni

re_clustering_draft <- function(features_by_cells, n_comp = 10, file1= "ElbowPlot.jpeg", file2 = "clustering.jpeg"){
	
	features_by_cells <- CreateSeuratObject(counts = features_by_cells, min.cells = 0, min.features = 0)
	all.genes <- rownames(features_by_cells)
	
	features_by_cells <- ScaleData(features_by_cells, features = all.genes)
	features_by_cells <- RunPCA(features_by_cells, features = all.genes)
	
	jpeg(file1, width=180, height=180, res=300, units="mm")
	ElbowPlot(features_by_cells)
	dev.off()
	
	features_by_cells <- FindNeighbors(features_by_cells, dims = 1:n_comp)
	features_by_cells <- FindClusters(features_by_cells)
	features_by_cells <- RunUMAP(features_by_cells, dims = 1:n_comp)
	
	jpeg(file2, width=180, height=180, res=300, units="mm")
	DimPlot(features_by_cells, reduction = "umap")
	dev.off()
	
	return(features_by_cells)
	
}
