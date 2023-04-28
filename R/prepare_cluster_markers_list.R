#' Prepare cluster markers list
#' Prepare cluster markers list for inter_dataset_comparison
#' @param marker_table_list list of marker tables, as produced by Seurat::FindAllMarkers
#' @param effect_size name of effect size column
#' @param effect_size_thr threshold over effect size
#' @param sig_column column with significance score
#' @param sig_thr threshold over effect significance score
#' @param id_column column with IDs
#' @param cluster column with cell cluster
#' @param top maximum number of markers per cluster
#' @export

prepare_cluster_markers_list <- function(marker_table_list=NULL, effect_size="avg_logFC", effect_size_thr=0.25, sig_column="p_val_adj", sig_thr=0.1, id_column=0, cluster="cluster", top=500){
  
  
  for(i in 1:length(marker_table_list)){
    
    #marker selection
    marker_table_list[[i]] <- marker_table_list[[i]][(marker_table_list[[i]][, sig_column] < sig_thr) & marker_table_list[[i]][, effect_size] > effect_size_thr, ]
    
    if(nrow(marker_table_list[[i]])>0){
      
      #split in cluster
      if(id_column==0){
        marker_table_list[[i]] <- split(rownames(marker_table_list[[i]]), marker_table_list[[i]][, cluster])  
      }else{
        marker_table_list[[i]] <- split(marker_table_list[[i]][, id_column], marker_table_list[[i]][, cluster])  
      }
      
      #limit each cluster to top elements
      marker_table_list[[i]] <- lapply(marker_table_list[[i]], function(i_list) lapply(i_list, function(i_vector) i_vector[1:min(top, length(i_vector))]))
      
      #create unique dataset markers names
      names(marker_table_list[[i]]) <- paste0(names(marker_table_list)[i], "_M", names(marker_table_list[[i]]))
      
    }
    
  }
  
  #from a list of lists to a gene set list
  gsl <- marker_table_list[[1]]
  for(i in 2:length(marker_table_list)){
    gsl <- c(gsl, marker_table_list[[i]])
  }
  marker_table_list <- gsl
  
  return(gsl)
  
}

