#' plot_umap_colored_features
#'
#' Generate a UMAP visualization colored by feature distributions 
#' @param Seu_obj seurat object, object with saved dimension reduction components calculate on genes by cells matrix
#' @param scMuffinList scMuffinList object
#' @param feature_id feature id
#' @param scale_feature whether to scale features or not
#' @param feature_breaks breaks for feature value distribution
#' @param ... further arguments to plot_umap
#' @param adj_outliers logical, whether to adjust the feature scores, removing outliers
#' @param min_cells min number of cells in which the feature must have a non-zero value
#' @param n_intervals number of intervals for the palette
#' @import pals
#' @description Generate a UMAP visualization colored by feature distributions
#' @export

plot_umap_colored_features <- function(Seu_obj=NULL, scMuffinList=NULL, feature_id=NULL, scale_feature=TRUE, n_intervals=10, adj_outliers=FALSE, min_cells=10, ...){
	
	# if(!dir.exists(dir)){
	# 	dir.create(dir)
	# }
	
	#feature_data <- features$df
	feature_data <- scMuffinList[[feature_id]]$summary
	feature_data <- feature_data[match(colnames(Seu_obj), rownames(feature_data)), ]
	
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
						
						#cut the feature id separately for <=0 and >0 in a total of 10 intervals
						#reorder cells according to feature data
						idx_neg <- feature_data_i <= 0
						idx_pos <- feature_data_i > 0
						md <- c(setNames(ggplot2::cut_number(feature_data_i[idx_neg], n_intervals/2), colnames(feature_data)[idx_neg]), setNames(ggplot2::cut_number(feature_data_i[feature_data_i>0], n_intervals/2), colnames(feature_data)[idx_pos]))
						md <- md[match(colnames(feature_data), names(md))]
						
						#add the metadata
						Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=md, col.name=colnames(feature_data)[i])
						if(any(is.na(Seu_obj@meta.data))){
							print(i)
							print(colnames(feature_data)[i])
							print(summary(Seu_obj@meta.data))
						}

						#print(table(Seu_obj@meta.data[, ncol(Seu_obj@meta.data)]))
						
						#pal <- rev(brewer.pal(11, "RdYlBu"))
						pal <- rev(pals::brewer.rdylbu(n_intervals))
						plot_umap(Seu_obj, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], cols=pal, ...)
						
					}else{
						
					  #if only positive values
						Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=feature_data_i[match(rownames(feature_data), colnames(Seu_obj))], col.name=colnames(feature_data)[i])
						plot_umap(Seu_obj, file = paste0(dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg"), group.by = colnames(feature_data)[i], feature_plot=TRUE, ...)
						
					}
				  
				}else{
				  
				  #if not numeric...
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
