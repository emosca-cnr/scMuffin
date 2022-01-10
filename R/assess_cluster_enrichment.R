#' assess_cluster_enrichment
#' @param features feature list
#' @param clusterings clustering list
#' @param meta_clusters boolean, TRUE if the assessment involves meta clusters
#' @param write_output boolean, TRUE to write the output in the output dir
#' @param out_dir output dir
#' @export
#' @importFrom utils write.table

assess_cluster_enrichment <- function(features, clusterings, meta_clusters=FALSE, write_output=TRUE, out_dir="./"){
  
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
      ans$cluster_gsea_res <- cluster_gsea(X_num, setNames(clusterings$meta_cl, clusterings$cell_id))
      
    }
    
    ##### FACTOR -> Fisher exact test
    if(ncol(X_fact)>0){
      cat("ORA for ", ncol(X_fact), "features\n")
      ans$cluster_hyper_res <- cluster_hyper(X_fact, setNames(clusterings$meta_cl, clusterings$cell_id))
      
    }
    
  }else{
    
    ##### NUMERCIC -> GSEA
    if(ncol(X_num)>0){
      cat("GSEA for ", ncol(X_num), "features\n")
      
      ans$cluster_gsea_res <- vector("list", ncol(clusterings))
      names(ans$cluster_gsea_res) <- colnames(clusterings)
      
      for(i in 1:ncol(clusterings)){
        ans$cluster_gsea_res[[i]] <- cluster_gsea(X_num, setNames(clusterings[, i], rownames(clusterings)))
      }
      
      
    }
    
    
  }
  
  ##### FACTOR -> Fisher exact test
  ##### FACTOR -> Fisher exact test
  if(ncol(X_fact)>0){
    cat("ORA for ", ncol(X_fact), "features\n")
    
    ans$cluster_hyper_res <- vector("list", ncol(clusterings))
    names(ans$cluster_hyper_res) <- colnames(clusterings)

    for(i in 1:ncol(clusterings)){
      ans$cluster_hyper_res[[i]] <- cluster_hyper(X_fact, setNames(clusterings[, i], rownames(clusterings)))
    }
    
    
    
  }
  
  
  
  return(ans)
  
}