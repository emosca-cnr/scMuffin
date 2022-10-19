#' plot_umap_colored_features
#'
#' Generate a UMAP visualization colored by feature distributions 
#' @param Seu_obj seurat object, object with saved dimension reduction components calculate on genes by cells matrix
#' @param features matrix, features by cells matrix
#' @param dir output directory
#' @param scale_feature whether to scale features or not
#' @param feature_breaks breaks for feature value distribution
#' @param ... further arguments to plot_umap
#' @param adj_outliers logical, whether to adjust the feature scores, removing outliers
#' @param min_cells min number of cells in which the feature must have a non-zero value
#' @param neg_eps small negative value to define the color for null values around zero, default to -0.1
#' @param pos_eps small positive value to define the color for null values around zero, default to =0.1
#' @import RColorBrewer
#' @description Generate a UMAP visualization colored by feature distributions
#' @export

plot_umap_colored_features <- function(Seu_obj, features, dir="./", scale_feature=TRUE, feature_breaks=NULL, adj_outliers=FALSE, neg_eps=-0.1, pos_eps=0.1, min_cells=10, ...){
	
	if(!dir.exists(dir)){
		dir.create(dir)
	}
	
	feature_data <- features$df
	#feature_data[is.na(feature_data)] <- 0
	#feature_data <- feature_data[, colSums(abs(feature_data))>0, drop=F]
	
	if(ncol(feature_data) > 0){
		
		
		for(i in 1:ncol(feature_data)){
			
			feature_data_i <- feature_data[, i]
			
			if(sum(!is.na(feature_data_i))>min_cells){
				
				if(is.numeric(feature_data_i)){
					
					if(adj_outliers){
						feature_data_i <- adj_outliers_col(feature_data_i)
					}
					
					if(scale_feature){
						feature_data_i <- scale(feature_data_i)
					}
					
					if(any(feature_data_i < 0, na.rm = T)){
						
						feature_data_i[is.na(feature_data_i)] <- 0
						
						if(is.null(feature_breaks)){
							fb <- c(seq(min(feature_data_i), neg_eps, length.out = 5), 0, seq(pos_eps, max(feature_data_i), length.out = 6))
						}
						
						Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=cut(feature_data_i[match(rownames(feature_data), colnames(Seu_obj))], breaks = fb, include.lowest = TRUE, right = TRUE), col.name=colnames(feature_data)[i])
						if(any(is.na(Seu_obj@meta.data))){
							print(i)
							print(colnames(feature_data)[i])
							print(summary(Seu_obj@meta.data))
						}

						#print(table(Seu_obj@meta.data[, ncol(Seu_obj@meta.data)]))
						
						pal <- rev(brewer.pal(11, "RdYlBu"))
						plot_umap(Seu_obj, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], cols=pal, ...)
						
					}else{
						
						Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=feature_data_i[match(rownames(feature_data), colnames(Seu_obj))], col.name=colnames(feature_data)[i])
						plot_umap(Seu_obj, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], feature_plot=TRUE, ...)
						
					}
				}else{
					
					plot_umap(Seu_obj, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], ...)
					
				}
				
			}else{
				message("not enough values for ", colnames(feature_data)[i], "\n")
			}
			
			
		}
	}else{
		message("No available values in input data\n")
		
	}
	
}
