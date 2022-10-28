#' Calculate gene set score at cluster level
#' @description Calculate gene set score at cluster level
#' @param gs_scores_obj list of gene signature lists
#' @param partition_id cell cluster labels
#' @param ncells_min min number of cells required for the clculation of the average signature in the cluster
#' @param null_model TRUE to consider the permutations
#' @param alt alterative
#' @param test type of test
#' @return scMuffinList with updates in elements cluster_data and cluster_data_full
#' @export

calculate_gs_scores_in_clusters <- function(scMuffinList=NULL, partition_id=NULL, ncells_min = 5, null_model = TRUE, alt="g", test="t"){
  
  cat("Clusters...\n")
  print(table(setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))))
  
  res_signatures_clusters <- lapply(scMuffinList$gene_set_scoring$full, function(i_marker_res) gs_scores_in_clusters(i_marker_res, cell_clusters=setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions)), ncells_min = ncells_min, null_model = null_model, alt=alt, test=test))
  
  #signatures-by-clusters matrix
  SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
  
  SC_signatures_by_cluster_matrix <- SC_signatures_by_cluster_matrix[, match(levels(setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))), colnames(SC_signatures_by_cluster_matrix))]
  
  
  ### OUTPUT
  
  scMuffinList$cluster_data[[partition_id]] <- list(gene_set_scoring=
                                                      list(
                                                        summary=as.data.frame(t(SC_signatures_by_cluster_matrix)),
                                                        full=res_signatures_clusters
                                                      )
  )
  
  return(scMuffinList) 
  
}