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
#' @param ... arguments passed to calculate_signatures
#' @export
#' @import RColorBrewer
#' @importFrom circlize colorRamp2

quantify_samples_similarity <- function(gbc_1, gbc_2, clusters_1, clusters_2, cluster_markers_1, cluster_markers_2, genes_min=3, genes_max=500, mc.cores=2, null_model=TRUE, ncells_min=5, do_plot=TRUE, sample_name_1="S1", sample_name_2="S2", outfile="cluster_similarity.jpg", pal=NULL, cluster_rows = F, cluster_columns = F, ...){
  
  if(any(lengths(cluster_markers_1)) > genes_max | any(lengths(cluster_markers_1) < genes_min)){
    message("number of genes beyond the chosen limits in at least one signature\n")
  }
  
  clusters_1 <- clusters_1[match(colnames(gbc_1), names(clusters_1))]
  clusters_2 <- clusters_2[match(colnames(gbc_2), names(clusters_2))]
  
  cl_names <- names(cluster_markers_1)
  cl_names_2 <- names(cluster_markers_2)
  
  
  signatures <- prepare_signatures(custom_signatures = list(cm=cluster_markers_1), genes = rownames(gbc_1), genes_min = genes_min, genes_max = genes_max)
  print(lengths(signatures$cm))
  
  
  res_sig <- mclapply(signatures, function(x) calculate_signatures(gbc_2, signatures = x, mc.cores=mc.cores/2, ...), mc.cores = mc.cores/2)
  
  #signature cluster median
  res_sig_cl <- mclapply(res_sig, function(x) calculate_signatures_clusters(sign_list = x$full, cell_clusters = clusters_2, null_model = null_model, ncells_min = ncells_min))
  
  
  print(lengths(cluster_markers_2))
  if(any(lengths(cluster_markers_2))>genes_max | any(lengths(cluster_markers_2) < genes_min)){
    message("number of genes beyond the chosen limits in at least one signature\n")
  }
  signatures_2 <- prepare_signatures(custom_signatures = list(cm=signatures_2), genes = rownames(gbc_2), genes_min = genes_min)
  
  res_sig_2 <- mclapply(signatures_2, function(x) calculate_signatures(gbc_1, signatures = x, mc.cores=mc.cores/2, ...), mc.cores = mc.cores/2)
  
  #signature cluster median
  res_sig_cl_2 <- mclapply(res_sig_2, function(x) calculate_signatures_clusters(sign_list = x$full, cell_clusters = clusters_1, null_model = null_model, ncells_min = ncells_min))
  
  
  #check for possible missing clusters
  m1 <- res_sig_cl$cm$signatures_by_clusters
  if(any(!cl_names %in% rownames(m1))){
    missing_rows <- as.character(cl_names[!cl_names %in% rownames(m1)])
    m1 <- rbind(m1, matrix(0, nrow = length(missing_rows), ncol = ncol(m1), dimnames = list(missing_rows, colnames(m1))))
  }
  m1 <- m1[order(as.numeric(rownames(m1))), order(as.numeric(colnames(m1)))]
  
  m2 <- res_sig_cl_2$cm$signatures_by_clusters
  if(any(!cl_names_2 %in% rownames(m2))){
    missing_rows <- as.character(cl_names_2[!cl_names_2 %in% rownames(m2)])
    m2 <- rbind(m2, matrix(0, nrow = length(missing_rows), ncol = ncol(m2), dimnames = list(missing_rows, colnames(m2))))
  }
  m2 <- m2[order(as.numeric(rownames(m2))), order(as.numeric(colnames(m2)))]
  
  print(m1)
  print(m2)
  
  clust_sim <- (m1 + t(m2)) / 2
  
  
  if(do_plot){
    
    if(is.null(pal)){
      pal <- RColorBrewer::brewer.pal(11, "PuOr")
    }
    
    extremes <- boxplot.stats(as.numeric(clust_sim))$stats
    extremes <- max(abs(extremes))
    col_fun <- circlize::colorRamp2(c(-extremes, 0, extremes), c(pal[1], pal[6], pal[11]))
    
    jpeg(outfile, width = 200, height = 180, res=300, units="mm")
    hm <- Heatmap(clust_sim, cluster_rows = cluster_rows, cluster_columns = cluster_columns, col = col_fun, name = "Similarity", rect_gp = gpar(col = "white", lwd = 2), row_title = sample_name_1, column_title = sample_name_2)
    draw(hm)
    dev.off()
    
  }
  
  return(list(clust_sim=clust_sim, markers_1=signatures$cm, markers_2=signatures_2$cm))
  
}