#' Cluster by features
#' 
#' Clustering of the features_by_cell matrix
#' @param features_by_cells matrix, features by cells matrix 
#' @param n_comp numeric, Dimensions of reduction to use as input
#' @param cnv TRUE/FALSE set it to TRUE for clustering CNV results 
#' @return features_by_cells Seurat Object, object with saved dimension reduction components calculate on features by cells matrix
#' @import Seurat
#' @export

cluster_by_features <- function(features, n_comp = 10, cnv=FALSE, plot_umap=FALSE, out_dir="./", scale_features=TRUE, ...){
	
	if(!dir.exists(out_dir)){
		dir.create(out_dir)
	}
	
	if(cnv){
		features_by_cells <- features
	}else{
		features_by_cells <- t(features$df)
		features_by_cells[is.na(features_by_cells)] <- 0
	}
	
	n_comp <- min(n_comp, nrow(features_by_cells))
	
	features_by_cells <- CreateSeuratObject(counts = features_by_cells, min.cells = 0, min.features = 0)
	all.genes <- rownames(features_by_cells)
	#scale data
	if(scale_features){
		features_by_cells <- ScaleData(features_by_cells, features = all.genes)
	}
	
	features_by_cells <- RunPCA(features_by_cells, features = all.genes)
	
	features_by_cells <- FindNeighbors(features_by_cells, dims = 1:n_comp)
	features_by_cells <- FindClusters(features_by_cells)
	features_by_cells <- RunUMAP(features_by_cells, dims = 1:n_comp)
	
	if(plot_umap){
		grDevices::jpeg(paste0(out_dir, "/feature_umap.jpg"), width = 180, height = 180, res=300, units="mm")
		res <- DimPlot(features_by_cells, ..., combine = F)
		plot(res[[1]])
		dev.off()
	}
	if(cnv){
		features_by_cells <- BuildClusterTree(features_by_cells, features = all.genes, reorder = T)
		hc_cells <- as.hclust(features_by_cells@tools$BuildClusterTree)
		ans <- list(clusters=features_by_cells@active.ident, hc=hc_cells, sobj=features_by_cells)
	}else{
		ans <- list(clusters=features_by_cells@active.ident, sobj=features_by_cells)
	}

	return(ans)
	
}
