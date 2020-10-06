#' scMuffin: main function that runs the whole pipeline
#' @param genes_by_cells seurat_object_data, genes-by-cells input matrix
#' @param custom_signatures list of gene sets
#' @import parallel

scMuffin <- function(genes_by_cells, cell_clustering=NULL, custom_signatures=NULL){


  #Seurat clustering @Noemi


  #signature expression @Ettore
  if(!is.null(custom_signatures)){
    signatures <- c(signatures, custom_signatures)
  }
  print(signatures)

  data_bins <- sc_data_bin(genes_by_cells, nbins = 25, use.log = TRUE, mc.cores=2)

  #parallel
  res_signatures <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score_in_clusters(i_marker_set, seurat_object_data = genes_by_cells, seurat_object_ident = seurat_clusters, perc=0.25, bins = data_bins, k=100, alt = "two.sided", test="t", nmark_min = 5, ncells_min=10), mc.cores = mc.cores)

  #Signalling entropy rate (SR) @Noemi
  #Potency states (LandScent): labels @Noemi
  #Diffusion pseudotime (DPT) (Destiny): @Noemi
  output_landscent <- landscent_sr(as.matrix(genes_by_cells@assays$RNA@data), mc.cores=mc.cores)

  #Paride score @Noemi
  exp_rate_score <- exp_rate(as.matrix(genes_by_cells@assays$RNA@counts))
  
  #Cell cycle state TO DO

  #CNV @Valentina
  cnv_res <- calculate_CNV(as.matrix(genes_by_cells@assays$RNA@data)[, 1:500], mc.cores = mc.cores)
  

  #merge everithing


  #re-clustering


  #Plots and Tables as outputs


  return(res)

}
