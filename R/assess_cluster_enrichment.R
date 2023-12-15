#' Assess enrichment with ORA or CSEA
#' @param scMuffinList scMuffinList list
#' @param feature_name the names of the feature that should be considered. It must be one of names(scMuffinList)
#' @param partition_id identifier of the partition to be considered
#' @param min.cells.feature minimum number of cells in which a feature must have value
#' @param min.cells.cluster minimum number of cells of a cluster
#' @param mc.cores number of cores
#' @param csea.k number of CSEA permutations
#' @param min.k minimum number of valid permutations to support empirical nulls
#' @param fract_min only cluster of size less or equal to this fraction of cell with not null feature values will be analysed
#' @description Assess cluster enrichment using ORA for categorical features and CSEA for numeric features.
#' @details
#' The output of CSEA is a table with statistics for every tested gene set. 
#' The output of ORA is composed of a series of tables with enrichment results, one for every possible categorical value. See [extract_cluster_enrichment_table] to extract summary table from CSEA and ORA results.
#'
#' @return scMuffinList with CSEA or ORA elements under scMuffinList$cluster_data for the considered partition
#' @importFrom stats setNames
#' @export

assess_cluster_enrichment <- function(scMuffinList=NULL, feature_name=NULL, partition_id=NULL, min.cells.feature=100, min.cells.cluster=10, fract_min=0.2, mc.cores=1, csea.k=99, min.k=10){
  
  if(length(scMuffinList[[feature_name]]$summary) == 0){
    stop("Can't find scMuffinList[[feature_name]]$summary\n")
  }
  if(!any(colnames(scMuffinList$partitions) == partition_id)){
    stop("Can't find any parition named ", partition_id, "\n")
  }
  
  
  X <- scMuffinList[[feature_name]]$summary
  X_clusters <- setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))
  
  ans <- vector("list", 2)
  names(ans) <- c("cluster_csea_res", "cluster_hyper_res")
  

  ##### NUMERCIC -> CSEA
  if(is.numeric(X[, 1])){
    
    cat("CSEA for ", ncol(X), "features\n")
    
    ans$cluster_csea_res <- cluster_csea(X, cell_clusters = X_clusters, min.cells.feature=min.cells.feature, min.cells.cluster=min.cells.cluster, mc.cores=mc.cores, csea.k=csea.k, min.k=min.k, fract_min=fract_min)
    
    if(length(scMuffinList$cluster_data[partition_id])==0){ # nor ORA neither CSEA
      
      scMuffinList$cluster_data[[partition_id]] <- list(CSEA=ans$cluster_csea_res)
      
    }else{
      
      if(length(scMuffinList$cluster_data[partition_id]$CSEA)==0){ # ORA already there
        
        scMuffinList$cluster_data[[partition_id]]$CSEA <- ans$cluster_hyper_res
        
      }
      
     
      shared_names <- which(names(scMuffinList$cluster_data[[partition_id]]$CSEA) %in% names(ans$cluster_csea_res))
      
      if(length(shared_names)>0){
        scMuffinList$cluster_data[[partition_id]]$CSEA[shared_names] <- NULL
      }

      scMuffinList$cluster_data[[partition_id]]$CSEA <- c(scMuffinList$cluster_data[[partition_id]]$CSEA, ans$cluster_csea_res)

    }
    
  }
  
  
  ##### FACTOR -> Fisher exact test
  if(is.factor(X[, 1]) | is.character(X[, 1])){
    
    cat("ORA for ", ncol(X), "features\n")
    
    ans$cluster_hyper_res <- cluster_hyper(X, cell_clusters = X_clusters)
    
    if(length(scMuffinList$cluster_data[partition_id])==0){ # nor ORA neither CSEA
      
      scMuffinList$cluster_data[[partition_id]] <- list(ORA=ans$cluster_hyper_res)
      
    }else{
      
      if(length(scMuffinList$cluster_data[partition_id]$ORA)==0){ #CSEA already there
        
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