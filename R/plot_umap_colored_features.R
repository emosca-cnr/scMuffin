#' UMAP Visualization colored by features
#'
#' @description Generate UMAP visualizations for all the columns in the "summary" element of a feature available in scMuffinList. It requires a Seurat object with UMAP information.
#' @param Seu_obj Seurat object with UMAP data
#' @param scMuffinList scMuffinList object
#' @param feature_name feature name
#' @param scale_feature whether to scale feature values
#' @param adj_outliers logical, whether to adjust the feature scores, removing outliers
#' @param min_cells minimum number of cells in which the feature must have a non-zero value
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param image_format png or jpeg
#' @param out_dir output directory
#' @param ... further arguments to plot_umap
#' @importFrom pals brewer.rdylbu brewer.purples alphabet
#' @importFrom stats setNames
#' @export

plot_umap_colored_features <- function(Seu_obj=NULL, scMuffinList=NULL, feature_name=NULL, scale_feature=TRUE, adj_outliers=FALSE, min_cells=10, out_dir="./", width=180, height=180, units="mm", res=300, image_format="png", ...){
  
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  
  #feature_data <- features$df
  feature_data <- scMuffinList[[feature_name]]$summary
  feature_data <- feature_data[match(colnames(Seu_obj), rownames(feature_data)), , drop=FALSE]
  
  #feature_data[is.na(feature_data)] <- 0
  #feature_data <- feature_data[, colSums(abs(feature_data))>0, drop=F]
  
  if(ncol(feature_data) > 0){
    
    for(i in 1:ncol(feature_data)){
      
      out_file <- paste0(out_dir, "/umap_by_genes_col_feature_", colnames(feature_data)[i],".jpg")
      
      feature_data_i <- feature_data[, i]
      
      if(sum(!is.na(feature_data_i))>min_cells){
        
        if(is.numeric(feature_data_i)){
          
          if(adj_outliers){
            feature_data_i <- adj_outliers_col(feature_data_i)
          }
          
          if(scale_feature){
            feature_data_i <- scale(feature_data_i)
          }
          
          feature_data_i[is.na(feature_data_i)] <- 0
          
          if(any(feature_data_i < 0)){
            
            #cut the feature id separately for <=0 and >0 in a total of 10 intervals
            #reorder cells according to feature data
            idx_neg <- feature_data_i <= 0
            idx_pos <- feature_data_i > 0
            md <- c(setNames(ggplot2::cut_interval(feature_data_i[idx_neg], 5, dig.lab = 2), rownames(feature_data)[idx_neg]), setNames(ggplot2::cut_interval(feature_data_i[feature_data_i>0], 5, dig.lab = 2), rownames(feature_data)[idx_pos]))
            md <- md[match(rownames(feature_data), names(md))]
            md <- factor(md, levels=rev(levels(md)))
            
            #add the metadata
            Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=md, col.name=colnames(feature_data)[i])
            plot_umap(Seu_obj, file = out_file, group.by = colnames(feature_data)[i], cols=(pals::brewer.rdylbu(10)), width=width, height=height, units=units, res=res, image_format=image_format, ...)
            
          }else{
            
            #if only positive values
            md <- ggplot2::cut_interval(feature_data_i, 5, dig.lab = 2)
            md <- factor(md, levels=rev(levels(md)))
            Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=md, col.name=colnames(feature_data)[i])
            plot_umap(Seu_obj, file = out_file, group.by = colnames(feature_data)[i], cols=rev(pals::brewer.ylorrd(5)), width=width, height=height, units=units, res=res, image_format=image_format, ...)
            
          }
          
        }else{
          
          #if not numeric...
          n_colors <- length(levels(as.factor(feature_data_i)))
          pal <- setNames(pals::alphabet(n_colors), levels(as.factor(feature_data_i)))
          Seu_obj <- Seurat::AddMetaData(Seu_obj, metadata=feature_data_i[match(rownames(feature_data), colnames(Seu_obj))], col.name=colnames(feature_data)[i])
          plot_umap(Seu_obj, file = out_file, group.by = colnames(feature_data)[i], cols=pal, width=width, height=height, units=units, res=res, image_format=image_format, ...)
          
        }
        
      }else{
        message("not enough values for ", colnames(feature_data)[i], "\n")
      }
      
      
    }
  }else{
    message("No available values in input data\n")
    
  }
  
}
