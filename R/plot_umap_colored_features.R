#' plot_umap_colored_features
#'
#' Generate a UMAP visualization colored by feature distributions 
#' @param seurat_object seurat object, object with saved dimension reduction components calculate on genes by cells matrix
#' @param features matrix, features by cells matrix
#' @param n_colors numeric, number of colors to be used
#' @import RColorBrewer
#' @export

plot_umap_colored_features <- function(seurat_object, features, dir="./", scale_feature=TRUE, feature_breaks=NULL, ...){
	
	if(!dir.exists(dir)){
		dir.create(dir)
	}
	
	feature_data <- features$df
	feature_data[is.na(feature_data)] <- 0
	feature_data <- feature_data[, colSums(abs(feature_data))>0, drop=F]
	
	if(ncol(feature_data) > 0){
		
		if(scale_feature){
			feature_data <- apply(feature_data, 2, scale)
			rownames(feature_data) <- rownames(features$df)
		}
		
		for(i in 1:ncol(feature_data)){
			
			feature_data_i <- feature_data[, i]
			
			if(is.numeric(feature_data_i)){
				
				if(any(feature_data_i < 0)){
					
					if(is.null(feature_breaks)){
						fb <- c(seq(min(feature_data_i), -0.1, length.out = 5), 0, seq(0.1, max(feature_data_i), length.out = 6))
					}
					seurat_object <- Seurat::AddMetaData(seurat_object, metadata=cut(feature_data_i[match(rownames(feature_data), colnames(seurat_object))], breaks = fb, include.lowest = TRUE, right = TRUE), col.name=colnames(feature_data)[i])
					if(any(is.na(seurat_object@meta.data))){
						print(i)
						print(colnames(feature_data)[i])
						print(summary(seurat_object@meta.data))
					}
					pal <- rev(brewer.pal(11, "RdYlBu"))
					plot_umap(seurat_object, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], cols=pal, ...)
					
				}else{
					
					seurat_object <- Seurat::AddMetaData(seurat_object, metadata=feature_data_i[match(rownames(feature_data), colnames(seurat_object))], col.name=colnames(feature_data)[i])
					plot_umap(seurat_object, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], feature_plot=TRUE, ...)
					
				}
			}else{
				
				plot_umap(seurat_object, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], ...)
				
			}
			
			
		}
	}else{
		warning("No available values in input data\n")
		
	}
	
}
