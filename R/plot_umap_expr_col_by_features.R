#' plot_umap_expr_features
#'
#' Generate a UMAP visualization colored by feature distributions 
#' @param seurat_object seurat object, object with saved dimension reduction components calculate on genes by cells matrix
#' @param features_matrix matrix, features by cells matrix
#' @import RColorBrewer
#' @export

plot_umap_expr_features <- function(seurat_object, cells_by_features, dir="./"){
	
	for(i in 1:ncol(cells_by_features)){
		feature_data <- cells_by_features[match(colnames(seurat_object), rownames(cells_by_features)), i]
		
		if(is.numeric(feature_data)){
			
			if(any(feature_data < 0)){
				
				feature_breaks <- c(seq(min(feature_data), -0.1, length.out = 5), 0, seq(0.1, max(feature_data), length.out = 6))
				seurat_object@meta.data$feature <- cut(feature_data, breaks = feature_breaks)
				pal <- brewer.pal(11, "RdYlBu")
				
			}else{
				
				seurat_object@meta.data$feature <- cut(feature_data, 9)
				pal <- brewer.pal(9, "Purples")
				
			}
			
		}else{
			
			seurat_object@meta.data$feature <- as.factor(feature_data)
			pal <- rainbow(length(levels(seurat_object@meta.data$feature)))
			
		}
		
		#UMAP by genes, colored by feature
		plot_umap(seurat_object, paste0(dir, "/umap_by_expr_col_by_", colnames(cells_by_features)[i],".jpg"), color_by="feature", pal=pal)
		
	}
	
}
