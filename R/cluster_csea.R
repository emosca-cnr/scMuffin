#' Calculate cluster enrichment by means of CSEA
#' 
#' @param feature_values matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param min.cells.feature minimum number of cells in which a feature must have value
#' @param min.cells.cluster minimum number of cells of a cluster
#' @param mc.cores number of cores
#' @param csea.k number of permutations
#' @param min.k minimum number of valid permutations to support empirical nulls
#' @return list with two elements: gs_table and leading_edge. See [csea()]
#' @export
#' @description Calculate cluster enrichment by csea approach


cluster_csea <- function(feature_values=NULL, cell_clusters=NULL, min.cells.feature=100, min.cells.cluster=10, mc.cores=1, csea.k=99, min.k=10){
  
  
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
  csea_res <- csea(X, gsl, mc_cores_perm = mc.cores, ord.mode = rep(-1, ncol(X)), k = csea.k, min.size = min.cells.feature, min.k = min.k)
  
  # insert a row for excluded clusters
  if(length(excluded_clusters)>0){
    
     csea_res <- lapply(csea_res, function(x) rbind(x, data.frame(id=excluded_clusters, es=NA, p_val=NA, adj_p_val=NA, nes=NA, FDRq=NA, nperm=NA, tags=NA, tags_perc=NA, list_top=NA, list_top_perc=NA, lead_edge=NA, row.names=excluded_clusters, stringsAsFactors = F)))
    
  }

  return(csea_res)
  
}
