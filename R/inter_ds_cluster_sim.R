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

inter_ds_cluster_sim <- function(gbc_1, gbc_2, clusters_1, clusters_2, cluster_markers_1, cluster_markers_2, genes_min=3, genes_max=500, mc.cores=1, null_model=TRUE, ncells_min=5, do_plot=TRUE, dataset_name_1="S1", dataset_name_2="S2", outfile="cluster_similarity.jpg", pal=NULL, cluster_rows = FALSE, cluster_columns = FALSE, top_genes=FALSE, ...){
  
  if(any(lengths(cluster_markers_1)) > genes_max | any(lengths(cluster_markers_1) < genes_min)){
    message("number of genes beyond the chosen limits in at least one signature\n")
  }
  

  scMuffinList_1 <- create_scMuffinList(normalized=Matrix::as.matrix(gbc_1)) 
  scMuffinList_1 <- add_partitions(scMuffinList_1, clusters = clusters_1, partition_id = "cl")
  
  scMuffinList_2 <- create_scMuffinList(normalized=Matrix::as.matrix(gbc_2)) 
  scMuffinList_2 <- add_partitions(scMuffinList_2, clusters = clusters_2, partition_id = "cl")
  
  cl_names <- names(cluster_markers_1)
  cl_names_2 <- names(cluster_markers_2)
  
  #the markers of dataset 1 must be present in dataset 2
  signatures <- prepare_gsls(custom_gsls = list(cm=cluster_markers_1), genes = rownames(scMuffinList_2$genes_by_cells), genes_min = genes_min, genes_max = genes_max)
  cat("dataset 1\n")
  print(lengths(signatures$cm))

  if(is.numeric(top_genes)){
    message("Reducing signatures to ", top_genes, "\n")
    signatures$cm <- lapply(signatures$cm, function(x) x[1:min(length(x), top_genes)])
  }
  print(lengths(signatures$cm))
  
  
  #expression of markers of dataset 1 in dataset 2
  scMuffinList_2 <- calculate_gs_scores(scMuffinList_2, gs_list = signatures$cm, mc.cores=mc.cores, ...)

  #signature cluster median
  scMuffinList_2 <- calculate_gs_scores_in_clusters(scMuffinList = scMuffinList_2, partition_id = "cl", null_model = null_model, ncells_min = ncells_min)
  
  
  cat("dataset 2\n")
  print(lengths(cluster_markers_2))
  if(any(lengths(cluster_markers_2))>genes_max | any(lengths(cluster_markers_2) < genes_min)){
    message("number of genes beyond the chosen limits in at least one signature\n")
  }
  
  #markers of dataset 2 must be in dataset 1
  signatures_2 <- prepare_gsls(custom_gsls = list(cm=cluster_markers_2), genes = rownames(scMuffinList_1$genes_by_cells), genes_min = genes_min)
  
  if(is.numeric(top_genes)){
    message("Reducing signatures to ", top_genes, "\n")
    signatures_2$cm <- lapply(signatures_2$cm, function(x) x[1:min(length(x), top_genes)])
  }
  print(lengths(signatures_2$cm))
  
  
  #expression of markers of dataset 2 in dataset 1
  scMuffinList_1 <- calculate_gs_scores(scMuffinList_1, gs_list = signatures_2$cm, mc.cores=mc.cores, ...)
  
  #signature cluster median
  scMuffinList_1 <- calculate_gs_scores_in_clusters(scMuffinList = scMuffinList_1, partition_id = "cl", null_model = null_model, ncells_min = ncells_min)
  
  
  #check for possible missing clusters
  m1 <- as.matrix(scMuffinList_1$cluster_data$cl$gene_set_scoring$summary)
  if(any(!cl_names %in% rownames(m1))){
    missing_rows <- as.character(cl_names[!cl_names %in% rownames(m1)])
    m1 <- rbind(m1, matrix(0, nrow = length(missing_rows), ncol = ncol(m1), dimnames = list(missing_rows, colnames(m1))))
  }
  m1 <- m1[order(as.numeric(rownames(m1))), order(as.numeric(colnames(m1)))]
  
  m2 <- as.matrix(scMuffinList_2$cluster_data$cl$gene_set_scoring$summary)
  if(any(!cl_names_2 %in% rownames(m2))){
    missing_rows <- as.character(cl_names_2[!cl_names_2 %in% rownames(m2)])
    m2 <- rbind(m2, matrix(0, nrow = length(missing_rows), ncol = ncol(m2), dimnames = list(missing_rows, colnames(m2))))
  }
  m2 <- m2[match(colnames(m1), rownames(m2)), match(rownames(m1), colnames(m2))]
  
  print(m1)
  print(m2)
  
  clust_sim <- (m1 + t(m2)) / 2
  
  
  if(do_plot){
    
    if(is.null(pal)){
      pal <- pals::brewer.puor(11)[c(1, 6, 11)]
    }
    
    extremes <- boxplot.stats(as.numeric(clust_sim))$stats
    extremes <- max(abs(extremes))
    col_fun <- circlize::colorRamp2(c(-extremes, 0, extremes), pal)
    
    jpeg(outfile, width = 200, height = 180, res=300, units="mm")
    hm <- Heatmap(clust_sim, cluster_rows = cluster_rows, cluster_columns = cluster_columns, col = col_fun, name = "Similarity", rect_gp = gpar(col = "white", lwd = 2), row_title = dataset_name_1, column_title = dataset_name_2)
    draw(hm)
    dev.off()
    
  }
  
  return(list(clust_sim=clust_sim, m1=m1, m2=m2, markers_1=signatures$cm, markers_2=signatures_2$cm))
  
}