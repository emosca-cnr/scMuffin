#' plot_gene_cluster_colored_features
#'
#'
#' @export

plot_gene_cluster_colored_features <- function(seurat_object, features_by_cells, n_colors=11){
	
	for(i in 1:nrow(features_by_cells)){
		seurat_object@meta.data$feature <- cut(GetAssayData(features_by_cells, slot = "scale.data")[i, ], n_colors)

		plot_umap(seurat_object, paste0("umap_by_genes_col_feature_", rownames(features_by_cells)[i],".jpg"), color_by="feature")
	}
	
}
