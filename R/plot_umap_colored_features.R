#' plot_umap_colored_features
#'
#' Generate a UMAP visualization colored by feature distributions 
#' @param seurat_object seurat object, object with saved dimension reduction components calculate on genes by cells matrix
#' @param features_by_cells matrix, features by cells matrix
#' @param n_colors numeric, number of colors to be used
#' @import RColorBrewer
#' @export

plot_umap_colored_features <- function(seurat_object, features_by_cells, dir="./"){
	
	for(i in 1:nrow(features_by_cells)){
		feature_data <- GetAssayData(features_by_cells, slot = "scale.data")[i, ]
		
		if(is.numeric(feature_data) & any(feature_data < 0)){
			
			feature_breaks <- c(seq(min(feature_data), -0.1, length.out = 5), 0, seq(0.1, max(feature_data), length.out = 6))
			seurat_object@meta.data$feature <- cut(feature_data, breaks = feature_breaks)
			pal <- brewer.pal(11, "RdYlBu")
			
		}else{
			
			seurat_object@meta.data$feature <- cut(feature_data, 9)
			pal <- brewer.pal(9, "Purples")
			
		}
	
		#UMAP by genes, colored by feature
		plot_umap(seurat_object, paste0(dir, "/umap_by_genes_col_feature_", rownames(features_by_cells)[i],".jpg"), color_by="feature", pal=pal)
		
		#UMAP by feature, colored by feature
		features_by_cells@meta.data$feature <- seurat_object@meta.data$feature
		plot_umap(features_by_cells, paste0(dir, "/umap_by_feature_col_feature_", rownames(features_by_cells)[i],".jpg"), color_by="feature", pal=pal)
		
	}
	
}
