#' Calculate cluster enrichment by means of CSEA
#' 
#' @param feature_values matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param min.cells.feature minimum number of cells in which a feature must have value
#' @param min.cells.cluster minimum number of cells of a cluster
#' @param mc.cores number of cores
#' @param csea.k number of permutations
#' @param min.k.nes minimum number of not null NES values
#' @return list with two elements: gs_table and leading_edge. See [csea()]
#' @export
#' @description Calculate cluster enrichment by csea approach


cluster_csea <- function(feature_values=NULL, cell_clusters=NULL, min.cells.feature=100, min.cells.cluster=10, mc.cores=1, csea.k=99, min.k.nes=10){
  
  
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
  
  #CSEA process on features-by-cells
  gsl <- lapply(split(cell_cluster_ok, cell_cluster_ok), function(x) names(x))
  csea_res <- csea(X, gsl, mc_cores_perm = mc.cores, ord.mode = rep(-1, ncol(X)), k = csea.k, min.size = min.cells.feature, min.k.nes = min.k.nes)
  
  # insert a row for excluded clusters
  if(length(excluded_clusters)>0){
    
    #debug
    cat(">>>>>>DEBUG<<<<<<")
    
    #gs table
    #csea_res$gs_table <- rbind(csea_res$gs_table, data.frame(id=excluded_clusters, es=0, p_val=1, adj_p_val=1, nes=0, FDRq=1, row.names=excluded_clusters, stringsAsFactors = F)) #old
    csea_res$gs_table <- lapply(csea_res$gs_table, function(x) rbind(x, data.frame(id=excluded_clusters, es=0, p_val=1, adj_p_val=1, nes=0, FDRq=1, n_pos_perm=0, n_neg_perm=0, row.names=excluded_clusters, stringsAsFactors = F)))
    
    csea_res$gs_table <- lapply(csea_res$gs_table, function(x) x[match(names(cluster_size), rownames(x)), ])
    
    #leading_edge
    csea_res$leading_edge <- lapply(csea_res$leading_edge, function(x) rbind(x, data.frame(tags=0, tags_perc=0, list_top=0, list_top_perc=0, lead_edge=0, lead_edge_subset="", row.names=excluded_clusters, stringsAsFactors = F)))
                                    
    #csea_res$leading_edge <- csea_res$leading_edge[match(names(cluster_size), rownames(csea_res$leading_edge)), ] #old
    csea_res$leading_edge <- lapply(csea_res$leading_edge, function(x) x[match(names(cluster_size), rownames(x)), ])
    
  }

  return(csea_res)
  
}
