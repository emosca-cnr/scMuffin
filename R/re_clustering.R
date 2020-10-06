#' gene_set_score_in_clusters
#' @param matrix matrix of features
#' @param n_comp Dimensions of reduction to use as input
#' @return seurat_m Seurat Object
#' 
#' @import Seurat
#' @export
#' 
re_clustering <- function(matrix, n_comp = 10, file1= "ElbowPlot.jpeg", file2 = "clustering.jpeg"){
seurat_m <- CreateSeuratObject(counts = temp_new, min.cells = 0, min.features = 0)
all.genes <- rownames(seurat_m)

seurat_m <- ScaleData(seurat_m, features = all.genes)
seurat_m <- RunPCA(seurat_m, features = all.genes)
jpeg(file1, width=180, height=180, res=300, units="mm")
ElbowPlot(seurat_m)
dev.off()

seurat_m <- FindNeighbors(seurat_m, dims = 1:n_comp)
seurat_m <- FindClusters(seurat_m)
seurat_m <- RunUMAP(seurat_m, dims = 1:n_comp)

jpeg(file2, width=180, height=180, res=300, units="mm")
DimPlot(seurat_m, reduction = "umap")
dev.off()

return(seurat_m)
}
