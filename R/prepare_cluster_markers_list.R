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

prepare_cluster_markers_list <- function(marker_table_list=NULL, effect_size="avg_log2FC", effect_size_thr=0.25, sig_column="p_val_adj", sig_thr=0.1, id_column="gene", cluster="cluster", top=500){
  
  
  if(!all(unlist(lapply(marker_table_list, function(x) all(c(effect_size, sig_column, cluster) %in% colnames(x)))))){
    stop("Any of ", effect_size, " ", sig_column, " ", cluster, " was not found in any element of marker_table_list.\n")
  }
  
  if(id_column!=0){
    if(!all(unlist(lapply(marker_table_list, function(x) id_column %in% colnames(x))))){
      stop(id_column, " was not found in any element of marker_table_list.\n")
    }
  }
  
  #
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
      marker_table_list[[i]] <- lapply(marker_table_list[[i]], function(i_list) i_list[1:min(top, length(i_list))])
      
      #create unique dataset markers names
      names(marker_table_list[[i]]) <- paste0(names(marker_table_list)[i], "_M", names(marker_table_list[[i]]))
      
    }
    
  }
  
  #from a list of lists to a gene set list
  gsl <- marker_table_list[[1]]
  for(i in 2:length(marker_table_list)){
    gsl <- c(gsl, marker_table_list[[i]])
  }
  print(lengths(gsl))
  
  return(gsl)
  
}

