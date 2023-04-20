#' Calculate cluster enrichment by means of GSEA
#' 
#' @param feature_values matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param min.cells.feature minimum number of cells in which a feature must have value
#' @param min.cells.cluster minimum number of cells of a cluster
#' @param mc.cores number of cores
#' @param gsea.k number of permutations
#' @return list with two elements: gs_table and leading_edge. See [gsea()]
#' @export
#' @description Calculate cluster enrichment by gsea approach
#' @author Ettore Mosca

cluster_gsea <- function(feature_values=NULL, cell_clusters=NULL, min.cells.feature=100, min.cells.cluster=10, mc.cores=1, gsea.k=99, min.k.nes=10){
  
  
  X <- as.matrix(feature_values)
  X[is.na(X)] <- 0
  
  #remove features with all null values
  idx_out <- colSums(X!=0) < min.cells.feature
  names_idx_out <- colnames(X)[idx_out]
  
  X <- X[, !idx_out, drop=F]
  if(nrow(X)<1){
    stop("not enough cells with values\n")
  }
  
  ## ensure that every cluster has at least min.cells.cluster with values !=0 
  cluster_size <- table(cell_clusters)
  cell_cluster_ok <- names(cluster_size)[cluster_size >= min.cells.cluster]
  excluded_clusters <- names(cluster_size)[cluster_size < min.cells.cluster]
  cell_cluster_ok <- cell_clusters[cell_clusters %in% cell_cluster_ok]
  
  #GSEA process on features-by-cells
  gsl <- lapply(split(cell_cluster_ok, cell_cluster_ok), function(x) names(x))
  gsea_res <- gsea(X, gsl, mc_cores_perm = mc.cores, ord.mode = rep(-1, ncol(X)), k = gsea.k, min.size = min.cells.feature, min.k.nes = min.k.nes)
  
  # insert a row for excluded clusters
  if(length(excluded_clusters)>0){
    
    #gs table
    gsea_res$gs_table <- rbind(gsea_res$gs_table, data.frame(id=excluded_clusters, es=0, p_val=1, adj_p_val=1, nes=0, FDRq=1, row.names=excluded_clusters, stringsAsFactors = F))
    gsea_res$gs_table <- gsea_res$gs_table[match(names(cluster_size), rownames(gsea_res$gs_table)), ]
    
    #leading_edge
    gsea_res$leading_edge <- rbind(gsea_res$leading_edge, data.frame(tags=0, tags_perc=0, list_top=0, list_top_perc=0, lead_edge=0, lead_edge_subset="", row.names=excluded_clusters, stringsAsFactors = F))
    gsea_res$leading_edge <- gsea_res$leading_edge[match(names(cluster_size), rownames(gsea_res$leading_edge)), ]
    
  }

  return(gsea_res)
  
}
