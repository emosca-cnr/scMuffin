#' CNV analysis
#' @param scMuffinList scMuffinList object
#' @param reference one-column data.frame with a reference expression profile; rownames must match those of genes_by_cell;
#' @param mc.cores number of cores;
#' @param wnd_size number of adjacent genes considered;
#' @param min_genes minimun number of genes expressed in a cell;
#' @param min_cells minimum numbver of cells in which a gene must be expressed;
#' @param expr_lim min and max values of relative expression; by default, lower and higher whiskers returned by grDevices::boxplot.stats will be used. Set to NA or anything other value such that length(expr_lim) != 2 to disable the removal of outliers.
#' @param center_cells whether to center cells. Default is TRUE
#' @param na.rm whether to remove 0 values in CNV estimation
#' @param center_genes whether to center genes or not. Default is FALSE. When using a reference it is useful to set this to TRUE
#' @param method "mean": subtract the average profile of the reference cluster to every cell (default); "min_max" (from Tirosh et al., DOI: 10.1126/science.aad0501) can be used but is not tested yet
#' @param z.score whether to use z-scores of the cluster median CNV profile, instead of the median itself.
#' @param eps absolute threshold to call CNV regions.
#' @param ... arguments passed to Seurat::FindClusters()
#' @param n_comp number of PCA components used for clustering
#' @param gene_ann optional data.frame with gene annotation with mandatory columns "Chromosome", "symbol" and "pos" (genomic location). If NULL gene locations will be collected from org.Hs.eg.db
#' @export 

CNV_analysis <- function(scMuffinList=NULL, mc.cores=1, reference = NULL, min_cells = 100, min_genes = 200, wnd_size = 100, center_cells = TRUE, center_genes = FALSE, expr_lim = FALSE, method="mean", na.rm=FALSE, z.score=FALSE, eps=NULL, n_comp=10, gene_ann=NULL, ...){
  
  cat("IMPORTANT: CNV analysis requires gene symbols as gene identifiers!\n")
  
  cnv_res <- calculate_CNV(scMuffinList$normalized, mc.cores = mc.cores, reference = reference, min_cells = min_cells, min_genes = min_genes, wnd_size = wnd_size, center_cells = center_cells, center_genes = center_genes, expr_lim = expr_lim, na.rm=na.rm, gene_ann=gene_ann)
  
  cnv_clustering <- cluster_by_features(cnv_res$CNV, n_comp=n_comp, ...) 

  if(!is.null(reference)){
    cnv_res$CNV <- apply_CNV_reference(cnv = cnv_res$CNV, cnv_clustering = cnv_clustering, reference="reference", method=method)
    cnv_res$ref_cluster <- cnv_res$CNV$ref_cluster
    cnv_res$CNV <- cnv_res$CNV$cnv
  }
  
  cnv_signal <- colSums(cnv_res$CNV^2)
  cnv_signal <- data.frame(CNV_score=as.numeric(cnv_signal), row.names = names(cnv_signal), stringsAsFactors = F)
  
  
  scMuffinList <- add_features(scMuffinList, name = "CNV", summary = cnv_signal, full=cnv_res)
  scMuffinList <- add_partitions(scMuffinList, partition_id = "CNV", clusters = cnv_clustering$clusters)

  scMuffinList$CNV$full$detected_cnv_regions <- detect_CNV_regions(scMuffinList=scMuffinList, z.score=z.score, eps=eps)
  
  return(scMuffinList)

}
