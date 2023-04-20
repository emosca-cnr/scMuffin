#' Inter-dataset cluster similarity
#' @param gbc_1 genes-by-cells expression matrix of dataset 1
#' @param gbc_2 genes-by-cells expression matrix of dataset 2
#' @param clusters_1 named vector with cell clustering fo dataset 1
#' @param clusters_2 named vector with cell clustering fo dataset 2
#' @param cluster_markers_1 list of markers for each cluster of dataset 1
#' @param cluster_markers_2 list of markers for each cluster of dataset 2
#' @param genes_min minimum number of genes required among the markers of a cluster
#' @param genes_max maximum number of genes required among the markers of a cluster
#' @param mc.cores number of cores
#' @param null_model whether to use or not the empirical null model. See calculate_signature
#' @param ncells_min minim number of cells in a cluster 
#' @param do_plot whether to plot or not the results
#' @param dataset_name_1 dataset 1 name
#' @param dataset_name_2 dataset 2 name
#' @param outfile out file name
#' @param pal palette for the output heatmap
#' @param cluster_rows whether to cluster or not the rows
#' @param cluster_columns whether to cluster or not the columns
#' @param top_genes If specified, only the first top_genes genes of every element of the lists cluster_markers_1 and cluster_markers_2 will be considered. This implies that the markers for each clusters are considered sorted by decreasing relevance. Default is FALSE, which
#' @param ... arguments passed to calculate_signatures
#' @export
#' @description Quantify the similarity between clusters of two datasets, on the basis of the average cluster marker expression
#' @import pals
#' @importFrom circlize colorRamp2
#' @return A list with:
#' \itemize{
#'   \item{clust_sim, matrix of clucster similarity;}
#'   \item{m1=m1}
#'   \item{m2=m2}
#'   \item{markers_1, markers of dataset 1;}
#'   \item{markers_2, markers of dataset 2;}
#' }

inter_dataset_comparison <- function(seu_obj_list=NULL, gsl=NULL, genes_min=3, genes_max=500, mc.cores=1, null_model=TRUE, ncells_min=5, do_plot=FALSE, outfile=NULL, pal=NULL, cluster_rows = FALSE, cluster_columns = FALSE, top_genes=FALSE,  ...){
  
  scMuffinList_list <- seu_obj_list
  score_matrix <- vector("list", length(scMuffinList_list))
  significance_matrix <- vector("list", length(scMuffinList_list))
  for(i in 1:length(seu_obj_list)){
    
    #scmuffin list  
    scMuffinList_list[[i]] <- create_scMuffinList(counts = GetAssayData(seu_obj_list[[i]], assay="RNA", slot="counts"), normalized=GetAssayData(seu_obj_list[[i]], assay="RNA", slot="data"))
    
    #partition
    scMuffinList_list[[i]] <- add_partitions(scMuffinList_list[[i]], clusters = setNames(paste0(names(seu_obj_list)[i], "_C", seu_obj_list[[i]]$seurat_clusters), names(seu_obj_list[[i]]$seurat_clusters)), partition_id = "C")
    
    #prepare gsl and calcualte gsl
    gsls <- prepare_gsls(custom_gsls = list(CM=gsl), scMuffinList = scMuffinList_list[[i]], genes_min = genes_min, genes_max = genes_max)
    scMuffinList_list[[i]] <- calculate_gs_scores(scMuffinList_list[[i]], gs_list = gsls$CM, nmark_min = genes_min, ncells_min = ncells_min)
    scMuffinList_list[[i]] <- calculate_gs_scores_in_clusters(scMuffinList_list[[i]], partition_id = "C", ncells_min = ncells_min)
    
    scMuffinList_list[[i]] <- assess_cluster_enrichment(scMuffinList_list[[i]], feature_name = "gene_set_scoring", partition_id = "C")
    
    score_matrix[[i]] <- scMuffinList_list[[i]]$cluster_data$C$gene_set_scoring$summary
    significance_matrix[[i]] <- extract_cluster_enrichment_table(scMuffinList_list[[i]],  partition_id = "C", type = "GSEA", quantity = "FDRq")
    
    #add missing cluster markers
    missing_marker_sets <- names(gsl)[!names(gsl) %in% colnames(score_matrix[[i]])]
    if(length(missing_marker_sets)>0){
      score_matrix[[i]] <- cbind(score_matrix[[i]], matrix(0, nrow = nrow(score_matrix[[i]]), ncol=length(missing_marker_sets), dimnames = list(rownames(score_matrix[[i]]), missing_marker_sets))) 
    }
    
    missing_marker_sets <- names(gsl)[!names(gsl) %in% colnames(significance_matrix[[i]])]
    if(length(missing_marker_sets)>0){
      significance_matrix[[i]] <- cbind(significance_matrix[[i]], matrix(1, nrow = nrow(significance_matrix[[i]]), ncol=length(missing_marker_sets), dimnames = list(rownames(significance_matrix[[i]]), missing_marker_sets))) 
    }
    
    score_matrix[[i]] <- score_matrix[[i]][, order(colnames(score_matrix[[i]]))]
    significance_matrix[[i]] <- significance_matrix[[i]][, order(colnames(significance_matrix[[i]]))]
    
    
  }
  rm(seu_obj_list, scMuffinList_list)
  
  score_matrix <- do.call(rbind, score_matrix)
  significance_matrix <- do.call(rbind, significance_matrix)
  
  zero_col_idx <- which(colSums(score_matrix)==0)
  if(length(zero_col_idx)>0){
    score_matrix <- score_matrix[, -zero_col_idx]
    significance_matrix <- significance_matrix[, -zero_col_idx]
  }
  
  ans <- list(score_matrix=score_matrix, significance_matrix=significance_matrix)
  
  
  return(ans)
  
}