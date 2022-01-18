#' extract_cluster_enrichment_table
#' @param clust_enrich_res cluster enrichment result
#' @param q_selection_criterion criteria for selection of tags from CSEA
#' @param q_selection_threshold threshold for selection of tags from CSEA
#' @param q_sort_crit criteria for sorting tags from CSEA
#' @param q_sort_desc sorting in decreasing order or not
#' @param only_pos_nes only positive nes?
#' @param c_selection_criterion criteria for selection of tags from ORA
#' @param c_selection_threshold threshold for selection of tags from ORA
#' @param c_sort_crit criteria for sorting tags from ORA
#' @param c_sort_desc sorting in decreasing order or ORA
#' @export

extract_cluster_enrichment_tags <- function(clust_enrich_res, q_selection_criterion="FDRq", q_selection_threshold=0.1, q_sort_crit="nes", q_sort_desc=FALSE, only_pos_nes=TRUE, c_selection_criterion="p_adj", c_selection_threshold=0.1, c_sort_crit="p", c_sort_desc=FALSE){
  
  
  ans <- vector("list", 2)
  names(ans) <- c("cluster_gsea_tags", "cluster_hyper_tags")
  
  en_table_sel <- extract_cluster_enrichment_table(clust_enrich_res, q_type=q_selection_criterion, c_type=c_selection_criterion)
  en_table_sort <- extract_cluster_enrichment_table(clust_enrich_res, q_type=q_sort_crit, c_type=c_sort_crit)
  
  if(only_pos_nes & q_sort_crit=="nes"){
    for(i in 1:length(en_table_sort$cluster_gsea_table)){
      idx_pos_nes <- en_table_sort$cluster_gsea_table[[i]] < 0
      en_table_sort$cluster_gsea_table[[i]][idx_pos_nes] <- 0 #negative nes ->0
      en_table_sel$cluster_gsea_table[[i]][idx_pos_nes] <- 1 #fdr equal to 1
    }
  }
  
  #selection of per-cluster tags
  en_table_sel$cluster_gsea_table <- lapply( en_table_sel$cluster_gsea_table, function(x) apply(x, 1, function(y) colnames(x)[y<q_selection_threshold]))
  
  #order
  for(i in 1:length(en_table_sort$cluster_gsea_table)){#clusterings
    n_clust <- length(en_table_sel$cluster_gsea_table[[i]])
    en_table_sort$cluster_gsea_table[[i]] <- split(t(apply(en_table_sort$cluster_gsea_table[[i]], 1, function(x) colnames(en_table_sort$cluster_gsea_table[[i]])[order(abs(x), decreasing = q_sort_desc)])), 1:n_clust)
  }
  
  #keep only those selected
  for(i in 1:length(en_table_sort$cluster_gsea_table)){#clusterings
    for(j in 1:length(en_table_sort$cluster_gsea_table[[i]])){#clusterings
      en_table_sort$cluster_gsea_table[[i]][[j]] <- en_table_sort$cluster_gsea_table[[i]][[j]][en_table_sort$cluster_gsea_table[[i]][[j]] %in% en_table_sel$cluster_gsea_table[[i]][[j]]]
    }
  }
  
  ans$cluster_gsea_tags <- en_table_sort$cluster_gsea_table
  
  #selection
  en_table_sel$cluster_hyper_table <- lapply( en_table_sel$cluster_hyper_table, function(x) lapply( x, function(y) apply(y, 1, function(z) colnames(y)[z<c_selection_threshold])))
  
  #order
  for(i in 1:length(en_table_sort$cluster_hyper_table)){#clusterings
    split_factor <- rownames(en_table_sort$cluster_hyper_table[[i]][[1]])
    en_table_sort$cluster_hyper_table[[i]] <- lapply(en_table_sort$cluster_hyper_table[[i]], function(x) split(t(apply(x, 1, function(y) colnames(x)[order(abs(y), decreasing = c_sort_desc)])), split_factor))
  }
  
  #keep only those selected
  for(i in 1:length(en_table_sort$cluster_hyper_table)){#clusterings
    for(j in 1:length(en_table_sort$cluster_hyper_table[[i]])){#features
      for(k in 1:length(en_table_sort$cluster_hyper_table[[i]][[j]])){#clusters
        en_table_sort$cluster_hyper_table[[i]][[j]][[k]] <- en_table_sort$cluster_hyper_table[[i]][[j]][[k]][en_table_sort$cluster_hyper_table[[i]][[j]][[k]] %in% en_table_sel$cluster_hyper_table[[i]][[j]][[k]] ]
      }
    }
  }
  #paste labels
  for(i in 1:length(en_table_sort$cluster_hyper_table)){#clusterings
    ans$cluster_hyper_tags[[i]] <- vector("list", length(en_table_sort$cluster_hyper_table[[i]][[1]])) #number of clustersfeatures
    names(ans$cluster_hyper_tags[[i]]) <- names(en_table_sort$cluster_hyper_table[[i]][[1]])
    for(j in 1:length(ans$cluster_hyper_tags[[i]])){
      ans$cluster_hyper_tags[[i]][[j]] <- lapply(en_table_sort$cluster_hyper_table[[i]], function(x) unlist(x[j]))
      ans$cluster_hyper_tags[[i]][[j]] <- paste(names(ans$cluster_hyper_tags[[i]][[j]]), lapply(ans$cluster_hyper_tags[[i]][[j]], function(x) paste0(x, collapse = "_")))
    }
  }
  
  return(ans)
  
}
