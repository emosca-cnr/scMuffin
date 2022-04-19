#' assess_cluster_enrichment
#' @param features feature list
#' @param partitions clustering list
#' @param meta_clusters boolean, TRUE if the assessment involves meta clusters
#' @param write_output boolean, TRUE to write the output in the output dir
#' @param out_dir output dir
#' @export
#' @importFrom utils write.table

assess_cluster_enrichment <- function(features, partitions, meta_clusters=FALSE, write_output=TRUE, out_dir="./"){
  
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive=TRUE)
  }
  
  X_fact <- as.matrix(features$df[, features$type=="factor", drop=F])
  X_num <- as.matrix(features$df[, !features$type=="factor", drop=F])
  X_num[is.na(X_num)] <- 0
  
  ans <- vector("list", 2)
  names(ans) <- c("cluster_gsea_res", "cluster_hyper_res")
  
  if(meta_clusters){
    
    ##### NUMERCIC -> GSEA
    if(ncol(X_num)>0){
      
      cat("GSEA for ", ncol(X_num), "features\n")
      ans$cluster_gsea_res <- cluster_gsea(X_num, setNames(partitions$meta_cl, partitions$cell_id))
      
    }
    
    ##### FACTOR -> Fisher exact test
    if(ncol(X_fact)>0){
      cat("ORA for ", ncol(X_fact), "features\n")
      ans$cluster_hyper_res <- cluster_hyper(X_fact, setNames(partitions$meta_cl, partitions$cell_id))
      
    }
    
  }else{
    
    ##### NUMERCIC -> GSEA
    if(ncol(X_num)>0){
      cat("GSEA for ", ncol(X_num), "features\n")
      
      ans$cluster_gsea_res <- vector("list", ncol(partitions))
      names(ans$cluster_gsea_res) <- colnames(partitions)
      
      for(i in 1:ncol(partitions)){
        ans$cluster_gsea_res[[i]] <- cluster_gsea(X_num, setNames(partitions[, i], rownames(partitions)))
      }
      
      
    }
    
    
  }
  
  ##### FACTOR -> Fisher exact test
  ##### FACTOR -> Fisher exact test
  if(ncol(X_fact)>0){
    cat("ORA for ", ncol(X_fact), "features\n")
    
    ans$cluster_hyper_res <- vector("list", ncol(partitions))
    names(ans$cluster_hyper_res) <- colnames(partitions)

    for(i in 1:ncol(partitions)){
      ans$cluster_hyper_res[[i]] <- cluster_hyper(X_fact, setNames(partitions[, i], rownames(partitions)))
    }
    
    
    
  }
  
  
  
  return(ans)
  
}