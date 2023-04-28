#' CNV analysis
#' @param scMuffinList scMuffinList object
#' @param reference one-column data.frame with a reference expression profile; rownames must match those of genes_by_cell;
#' @param mc.cores number of cores;
#' @param wnd_size number of adjacent genes considered;
#' @param min_genes minimun number of genes expressed in a cell;
#' @param min_cells minimum numbver of cells in which a gene must be expressed;
#' @param expr_lim min and max values of relative expression; by default, lower and higher whiskers returned by grDevices::boxplot.stats will be used. Set to NA or anything other value such that length(expr_lim) != 2 to disbale the removal of outliers.
#' @param scale_cells whether to scale cells
#' @param na.rm whether to remove 0 values in CNV estimation
#' @param center_genes whether to center genes or not
#' @param method mean: subtract the average profile of the reference cluster to every cell; min_max see Tirosh et al.
#' @param z.score whether to use z-scores of the cluster median CNV profile, instead of the median itself.
#' @param eps absolute threshold to call CNV regions.
#' @export 

CNV_analysis <- function(scMuffinList=NULL, mc.cores=1, reference = NULL, min_cells = 100, min_genes = 100, wnd_size = 100, scale_cells = TRUE, center_genes = FALSE, expr_lim = FALSE, method="mean", na.rm=FALSE, z.score=FALSE, eps=NULL){
  
  cnv_res <- calculate_CNV(scMuffinList$normalized, mc.cores = mc.cores, reference = reference, min_cells = min_cells, min_genes = min_genes, wnd_size = wnd_size, scale_cells = scale_cells, center_genes = center_genes, expr_lim = expr_lim, na.rm=na.rm)
  
  cnv_clustering <- cluster_by_features(cnv_res$CNV, cnv=TRUE) 

  if(!is.null(reference)){
    cnv_res$CNV <- apply_CNV_reference(cnv = cnv_res$CNV, cnv_clustering = cnv_clustering, reference="reference", method=method)
    cnv_res$ref_cluster <- cnv_res$CNV$ref_cluster
    cnv_res$CNV <- cnv_res$CNV$cnv
  }
  
  cnv_signal <- colSums(cnv_res$CNV^2)
  cnv_signal <- data.frame(CNV_signal=as.numeric(cnv_signal), row.names = names(cnv_signal), stringsAsFactors = F)
  
  
  scMuffinList <- add_features(scMuffinList, name = "CNV", summary = cnv_signal, full=cnv_res)
  scMuffinList <- add_partitions(scMuffinList, partition_id = "CNV", clusters = cnv_clustering$clusters)

  scMuffinList$CNV$full$detected_cnv_regions <- detect_CNV_regions(scMuffinList=scMuffinList, z.score=z.score, eps=eps)
  
  return(scMuffinList)

}