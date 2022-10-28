#' Compare two samples on the basis of cluster markers
#' @param gbc_1 genes-by-cells expression matrix of sample 1
#' @param gbc_2 genes-by-cells expression matrix of sample 2
#' @param clusters_1 named vector with cell clustering fo sample 1
#' @param clusters_2 named vector with cell clustering fo sample 2
#' @param cluster_markers_1 list of markers for each cluster of sample 1
#' @param cluster_markers_2 list of markers for each cluster of sample 2
#' @param genes_min minimum number of genes required among the markers of a cluster
#' @param genes_max maximum number of genes required among the markers of a cluster
#' @param mc.cores number of cores
#' @param null_model whether to use or not the empirical null model. See calculate_signature
#' @param ncells_min minim number of cells in a cluster 
#' @param do_plot whether to plot or not the results
#' @param sample_name_1 sample 1 name
#' @param sample_name_2 sample 2 name
#' @param outfile out file name
#' @param pal palette for the output heatmap
#' @param cluster_rows whether to cluster or not the rows
#' @param cluster_columns whether to cluster or not the columns
#' @param top_genes If specified, only the first top_genes genes of every element of the lists cluster_markers_1 and cluster_markers_2 will be considered. This implies that the markers for each clusters are considered sorted by decreasing relevance. Default is FALSE, which
#' @param ... arguments passed to calculate_signatures
#' @export
#' @description Compare two samples on the basis of cluster markers
#' @import pals
#' @importFrom circlize colorRamp2

quantify_samples_similarity <- function(gbc_1, gbc_2, clusters_1, clusters_2, cluster_markers_1, cluster_markers_2, genes_min=3, genes_max=500, mc.cores=2, null_model=TRUE, ncells_min=5, do_plot=TRUE, sample_name_1="S1", sample_name_2="S2", outfile="cluster_similarity.jpg", pal=NULL, cluster_rows = FALSE, cluster_columns = FALSE, top_genes=FALSE, ...){
  
  if(any(lengths(cluster_markers_1)) > genes_max | any(lengths(cluster_markers_1) < genes_min)){
    message("number of genes beyond the chosen limits in at least one signature\n")
  }
  

  scMuffinList_1 <- list(genes_by_cells=Matrix::as.matrix(gbc_1))
  scMuffinList_1 <- add_partitions(scMuffinList_1, clusters = clusters_1, partition_id = "cl")
  
  scMuffinList_2 <- list(genes_by_cells=Matrix::as.matrix(gbc_2))
  scMuffinList_2 <- add_partitions(scMuffinList_2, clusters = clusters_2, partition_id = "cl")
  
  cl_names <- names(cluster_markers_1)
  cl_names_2 <- names(cluster_markers_2)
  
  #the markers of sample 1 must be present in sample 2
  signatures <- prepare_gsls(custom_gsls = list(cm=cluster_markers_1), genes = rownames(scMuffinList_2$genes_by_cells), genes_min = genes_min, genes_max = genes_max)
  cat("Sample 1\n")
  print(lengths(signatures$cm))

  if(is.numeric(top_genes)){
    message("Reducing signatures to ", top_genes, "\n")
    signatures$cm <- lapply(signatures$cm, function(x) x[1:min(length(x), top_genes)])
  }
  print(lengths(signatures$cm))
  
  
  #expression of markers of sample 1 in sample 2
  scMuffinList_2 <- calculate_gs_scores(scMuffinList_2, gs_list = signatures$cm, mc.cores=mc.cores, ...)

  #signature cluster median
  scMuffinList_2 <- calculate_gs_scores_in_clusters(scMuffinList = scMuffinList_2, partition_id = "cl", null_model = null_model, ncells_min = ncells_min)
  
  
  cat("Sample 2\n")
  print(lengths(cluster_markers_2))
  if(any(lengths(cluster_markers_2))>genes_max | any(lengths(cluster_markers_2) < genes_min)){
    message("number of genes beyond the chosen limits in at least one signature\n")
  }
  
  #markers of sample 2 must be in sample 1
  signatures_2 <- prepare_gsls(custom_gsls = list(cm=cluster_markers_2), genes = rownames(scMuffinList_1$genes_by_cells), genes_min = genes_min)
  
  if(is.numeric(top_genes)){
    message("Reducing signatures to ", top_genes, "\n")
    signatures_2$cm <- lapply(signatures_2$cm, function(x) x[1:min(length(x), top_genes)])
  }
  print(lengths(signatures_2$cm))
  
  
  #expression of markers of sample 2 in sample 1
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
      pal <- pals::brewer.puor(11)
    }
    
    extremes <- boxplot.stats(as.numeric(clust_sim))$stats
    extremes <- max(abs(extremes))
    col_fun <- circlize::colorRamp2(c(-extremes, 0, extremes), c(pal[1], pal[6], pal[11]))
    
    jpeg(outfile, width = 200, height = 180, res=300, units="mm")
    hm <- Heatmap(clust_sim, cluster_rows = cluster_rows, cluster_columns = cluster_columns, col = col_fun, name = "Similarity", rect_gp = gpar(col = "white", lwd = 2), row_title = sample_name_1, column_title = sample_name_2)
    draw(hm)
    dev.off()
    
  }
  
  return(list(clust_sim=clust_sim, m1=m1, m2=m2, markers_1=signatures$cm, markers_2=signatures_2$cm))
  
}