#' Calculate gene set score at cluster level
#' @description Calculate gene set score at cluster level
#' @param scMuffinList scMuffinList object
#' @param partition_id identifier of the partition to be used
#' @param ncells_min minimum number of cells required for the calculation of the average signature in the cluster
#' @param alt alterative passed to [wilcox.test()] or [t.test()]
#' @param test type of test: t to use [t.test()]; wrs to use [wilcox.test()]
#' @param fract_min only clusters with this fraction of cells with not null gene set score will be considered
#' @return scMuffinList with cluster level scores in `sMuffinList$cluster_data[[partition_id]]`. The element [summary] contains a clusters-by-gene sets table, while the element [full] the full result
#' @export

calculate_gs_scores_in_clusters <- function(scMuffinList=NULL, partition_id=NULL, ncells_min = 5, alt="g", test="t", fract_min=0.5){
  
  if(!any(colnames(scMuffinList$partitions) == partition_id)){
    stop("Can't find any parition named ", partition_id, "\n")
  }
  
  cat("Clusters...\n")
  print(table(setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))))
  
  res_signatures_clusters <- lapply(scMuffinList$gene_set_scoring$full, function(i_marker_res) gs_scores_in_clusters(i_marker_res, cell_clusters=setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions)), ncells_min = ncells_min, fract_min = fract_min, alt=alt, test=test))
  
  #signatures-by-clusters matrix
  SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
  
  #SC_signatures_by_cluster_matrix <- SC_signatures_by_cluster_matrix[, match(levels(setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))), colnames(SC_signatures_by_cluster_matrix))]
  
  cluster_levels <- levels(factor(setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))))
  SC_signatures_by_cluster_matrix <- SC_signatures_by_cluster_matrix[, match(cluster_levels, colnames(SC_signatures_by_cluster_matrix))]

  ### OUTPUT
  
  scMuffinList$cluster_data[[partition_id]]$gene_set_scoring <- list(
                                                        summary=as.data.frame(t(SC_signatures_by_cluster_matrix)),
                                                        full=res_signatures_clusters
                                                      )

  
  return(scMuffinList) 
  
}