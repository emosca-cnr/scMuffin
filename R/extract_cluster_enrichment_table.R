#' extract_cluster_enrichment_table
#' @param features feature list
#' @param clusterings clustering list
#' @param meta_clusters boolean, TRUE if the assessment involves meta clusters
#' @param write_output boolean, TRUE to write the output in the output dir
#' @param out_dir output dir
#' @export
#' @importFrom utils write.table
extract_cluster_enrichment_table <- function(clust_enrich_res, q_type="nes", c_type="er"){

  
  ans <- vector("list", 2)
  names(ans) <- c("cluster_gsea_table", "cluster_hyper_table")
  
  if(!is.null(q_type) & !is.null(clust_enrich_res$cluster_gsea_res)){
    ans$cluster_gsea_table <- clust_enrich_res$cluster_gsea_res
    for(i in 1:length(clust_enrich_res$cluster_gsea_res)){
      ans$cluster_gsea_table[[i]] <- do.call(cbind, lapply(clust_enrich_res$cluster_gsea_res[[i]]$gs_table, function(x) array(x[, colnames(x)==q_type], dimnames = list(x$id))))
    }
      
      
  }
  
  if(!is.null(c_type) & !is.null(clust_enrich_res$cluster_hyper_res)){
    ans$cluster_hyper_table <- clust_enrich_res$cluster_hyper_res
    
    for(i in 1:length(clust_enrich_res$cluster_hyper_res)){ #every item corresponds to an item of clusterings
      for(j in 1:length(clust_enrich_res$cluster_hyper_res[[i]])){ #every item corresponds to a feature
        ans$cluster_hyper_table[[i]][[j]] <- do.call(cbind, lapply(clust_enrich_res$cluster_hyper_res[[i]][[j]], function(x) setNames(x$clusters[, colnames(x$clusters)==c_type], x$clusters$id)))
      }
    }


  }

  return(ans)
  
}