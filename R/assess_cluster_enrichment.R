#' Assess enrichment with ORA or GSEA
#' @param scMuffinList scMuffinList list
#' @param feature_name the names of the feature that should be considered. It must be one of names(scMuffinList)
#' @param partition_id identifier of the partition to be considered
#' @description Assess cluster enrichment using ORA for categorical features and GSEA for numeric features.
#' @details
#' The output of GSEA is a table with statistics for every tested gene set (gs_table element) and a table with the results of the leading edge (leading_edge element). 
#' The output of ORA is composed of a series of tables with enrichment results, one for every possible categorical value. See [extract_cluster_enrichment_table] to extract summary table from GSEA and ORA results.
#'
#' @return scMuffinList with GSEA or ORA elements under scMuffinList$cluster_data for the considered partition
#' @export

assess_cluster_enrichment <- function(scMuffinList=NULL, feature_name=NULL, partition_id=NULL, meta_clusters=FALSE){
  
  
  X <- scMuffinList[[feature_name]]$summary
  X_clusters <- setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))
  
  ans <- vector("list", 2)
  names(ans) <- c("cluster_gsea_res", "cluster_hyper_res")
  
  # if(meta_clusters){
  #   
  #   ##### NUMERCIC -> GSEA
  #   if(ncol(X)>0){
  #     
  #     cat("GSEA for ", ncol(X), "features\n")
  #     ans$cluster_gsea_res <- cluster_gsea(X, setNames(partitions$meta_cl, partitions$cell_id))
  #     
  #   }
  #   
  #   ##### FACTOR -> Fisher exact test
  #   if(ncol(X)>0){
  #     cat("ORA for ", ncol(X), "features\n")
  #     ans$cluster_hyper_res <- cluster_hyper(X, cell_clusters = X_clusters)
  #     
  #   }
  #   
  # }else{
  
  ##### NUMERCIC -> GSEA
  if(is.numeric(X[, 1])){
    
    cat("GSEA for ", ncol(X), "features\n")
    
    ans$cluster_gsea_res <- cluster_gsea(X, cell_clusters = X_clusters)
    
    if(length(scMuffinList$cluster_data[partition_id])==0){ # nor ORA neither GSEA
      
      scMuffinList$cluster_data[[partition_id]] <- list(GSEA=ans$cluster_gsea_res)
      
    }else{
      
      if(length(scMuffinList$cluster_data[partition_id]$GSEA)==0){ # ORA already there
        
        scMuffinList$cluster_data[[partition_id]]$GSEA <- lans$cluster_hyper_res
        
      }
      
      
      shared_names <- which(names(scMuffinList$cluster_data[[partition_id]]$GSEA$gs_table) %in% names(ans$cluster_gsea_res$gs_table))
      if(length(shared_names)>0){
        scMuffinList$cluster_data[[partition_id]]$GSEA$gs_table[shared_names] <- NULL
        scMuffinList$cluster_data[[partition_id]]$GSEA$leading_edge[shared_names] <- NULL
      }
      
      scMuffinList$cluster_data[[partition_id]]$GSEA$gs_table <- c(scMuffinList$cluster_data[[partition_id]]$GSEA$gs_table, ans$cluster_gsea_res$gs_table)
      scMuffinList$cluster_data[[partition_id]]$GSEA$leading_edge <- c(scMuffinList$cluster_data_full[[partition_id]]$GSEA$leading_edge, ans$cluster_gsea_res$leading_edge) ##append 
    }
    
  }
  
  
  ##### FACTOR -> Fisher exact test
  if(is.factor(X[, 1])){
    
    cat("ORA for ", ncol(X), "features\n")
    
    ans$cluster_hyper_res <- cluster_hyper(X, cell_clusters = X_clusters)
    
    if(length(scMuffinList$cluster_data[partition_id])==0){ # nor ORA neither GSEA
      
      scMuffinList$cluster_data[[partition_id]] <- list(ORA=ans$cluster_hyper_res)
      
    }else{
      
      if(length(scMuffinList$cluster_data[partition_id]$ORA)==0){ #GSEA already there
        
        scMuffinList$cluster_data[[partition_id]]$ORA <- ans$cluster_hyper_res
        
      }else{
        
        #overwrite possible pre-existing elements with the same name
        shared_names <- which(names(scMuffinList$cluster_data[[partition_id]]$ORA) %in% names(ans$cluster_hyper_res))
        if(length(shared_names)>0){
          scMuffinList$cluster_data[[partition_id]]$ORA[shared_names] <- NULL
        }
        scMuffinList$cluster_data[[partition_id]]$ORA <- c(scMuffinList$cluster_data[[partition_id]]$ORA, ans$cluster_hyper_res)
        
        #ensure the same ordering by matching the rownames of all elments to those of the first element
        for(i in 2:length(scMuffinList$cluster_data[[partition_id]]$ORA)){
          for(j in 1:length(scMuffinList$cluster_data[[partition_id]]$ORA[[i]])){
            scMuffinList$cluster_data[[partition_id]]$ORA[[i]][[j]] <- scMuffinList$cluster_data[[partition_id]]$ORA[[i]][[j]][match(rownames(scMuffinList$cluster_data[[partition_id]]$ORA[[1]][[1]]), , rownames(scMuffinList$cluster_data[[partition_id]]$ORA[[i]][[j]])), ]
          }
        }
      }
    }
  }
  
  
  return(scMuffinList)
  
}