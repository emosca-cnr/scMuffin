#' Extract the top results of GSEA and ORA for each cluster
#' @param scMuffinList scMuffinList object
#' @param partition_id one among the partitions
#' @param GSEA_selection_criterion selection criteria for results of GSEA
#' @param GSEA_selection_threshold threshold for selection of tags from GSEA
#' @param ORA_selection_criterion selection criteria for results of ORA
#' @param ORA_selection_threshold threshold for selection of tags from ORA
#' @param n_max_per_cluster maximum number of tags per cluster
#' @description Extract the best results of GSEA and ORA for each cluster
#' @return Add or overwrite the "cluster_tags" list to `scMuffinList$cluster_data[[partition_id]]`
#' @export

extract_cluster_enrichment_tags <- function(scMuffinList=NULL, partition_id=NULL, GSEA_selection_criterion="FDRq", GSEA_selection_threshold=0.25, only_pos_nes=TRUE, ORA_selection_criterion="p_adj", ORA_selection_threshold=0.1, n_max_per_cluster=3){
  
  clust_enrich_res <- scMuffinList$cluster_data[[partition_id]][c("GSEA", "ORA")]
  
  ans <- vector("list", 2)
  names(ans) <- c("cluster_gsea_tags", "cluster_hyper_tags")
  
  ##GSEA results
  GSEA_table <- extract_cluster_enrichment_table(scMuffinList = scMuffinList, partition_id = partition_id, type = "GSEA", quantity = GSEA_selection_criterion)
  if(only_pos_nes){
    GSEA_table_nes <- extract_cluster_enrichment_table(scMuffinList = scMuffinList, partition_id = partition_id, type = "GSEA", quantity = "nes")
    GSEA_table[GSEA_table_nes < 0] <- 1 #fdr equal to 1
  }
  
  #selection of per-cluster tags
  decreasing <- FALSE
  if(GSEA_selection_criterion == "nes"){
    decreasing <- TRUE
  }
  gsea_tags <- apply(GSEA_table, 1, function(i_row) intersect(colnames(GSEA_table)[order(i_row, decreasing = decreasing)], colnames(GSEA_table)[i_row < GSEA_selection_threshold]))
  gsea_tags <- lapply(gsea_tags, function(x) x[1:min(n_max_per_cluster, length(x))])
  
  #selection
  ORA_table <- extract_cluster_enrichment_table(scMuffinList = scMuffinList, partition_id = partition_id, type = "ORA", quantity = ORA_selection_criterion)
  ORA_table <- do.call(cbind, ORA_table)
  
  ora_tags <- apply(ORA_table, 1, function(i_row) intersect(colnames(ORA_table)[order(i_row, decreasing = decreasing)], colnames(ORA_table)[i_row < ORA_selection_threshold]))
  ora_tags <- lapply(ora_tags, function(x) x[1:min(n_max_per_cluster, length(x))])
  
  
  decreasing <- TRUE
  if(GSEA_selection_criterion %in% c("p", "p_adj")){
    decreasing <- FALSE
  }
  ora_tags <- apply(ORA_table, 1, function(y) intersect(colnames(ORA_table)[order(y, decreasing = decreasing)], colnames(ORA_table)[y<ORA_selection_threshold]))
  
  scMuffinList$cluster_data[[partition_id]]$cluster_tags <- list(GSEA=gsea_tags, ORA=ora_tags)
  
  return(scMuffinList)
  
}
