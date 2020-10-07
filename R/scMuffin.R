#' scMuffin: main function that runs the whole pipeline
#' @param genes_by_cells seurat_object_data, genes-by-cells input matrix
#' @param custom_signatures list of gene sets
#' @import parallel

scMuffin <- function(genes_by_cells, custom_signatures=NULL, mc.cores=2){


	##################	SIGNATURES #################
	cat("Calcuting gene signature scores...\n")
	
	if(!is.null(custom_signatures)){
    signatures <- c(signatures, custom_signatures)
  }
  cat("# of signatures: ", length(signatures), "\n")

  #dataset bins
  data_bins <- sc_data_bin(as.matrix(genes_by_cells@assays$RNA@data), nbins = 25, use.log = TRUE)

  #signatures-by-cell matrix
  res_signatures <- mclapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=100, nmark_min = 5, ncells_min = 5), mc.cores = mc.cores)

  #gene-set score per cluster list
  res_signatures_clusters <- lapply(res, function(i_marker_res) gene_set_score_in_clusters(i_marker_res$score_table, genes_by_cells@active.ident, ncells_min = 5))
  
  #signatures-by-clusters matrix
  SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))

  #output
  heatmap_signatures(SC_signatures_by_cluster_matrix)
  write.table(SC_signatures_by_cluster_matrix, file="signatures_by_clusters.txt", sep = "\t", row.names = T, col.names = NA)
  
  ##################	Signalling entropy rate (SR) @Noemi 	  ##################	
  ##################	Potency states (LandScent): labels @Noemi 	  ##################	
  ##################	Diffusion pseudotime (DPT) (Destiny): @Noemi   ##################	
  cat("Calcuting landscent-related scores...\n")
  output_landscent <- landscent_sr(as.matrix(genes_by_cells@assays$RNA@data), mc.cores=mc.cores)

  
  ##################	EXPRESSION RATE   ##################	
  exp_rate_score <- exp_rate(as.matrix(genes_by_cells@assays$RNA@counts))
  
  ##################	Cell cycle state   ##################	

  
  ##################	CNV @Valentina   ##################	
  cnv_res <- calculate_CNV(as.matrix(genes_by_cells@assays$RNA@data)[, 1:500], mc.cores = mc.cores)
  
  ngenes_chrom <- unlist(lapply(cnv_res, nrow)) # number of genes per chromosome
  cnv_res <- preprocess_for_heatmap_CNV(cnv_res)
  heatmap_CNV_clusters <- heatmap_CNV(cnv_res, ngenes_chrom)

  ##################	merge everithing   ##################	
  features_by_cells <- merge_matrix(signatures_by_cells = res_signatures, expr_score = exp_rate_score, cnv=cnv_res, output_landscent = output_landscent)

  feature_corr <- cor(features_by_cells, method="spearman") #will this work?
  
  
  ##################	re-clustering   ##################	
  features_by_cells <- re_clustering(features_by_cells)
  
  #### FINAL OUTPUTS ###

  #### Visual comparison between initial and final clusters ###
  plot_umap(genes_by_cells, "umap_genes.jpg")
  
  plot_umap(features_by_cells, "umap_features.jpg")
  
  #add initial clusters information in features_by_cells meta.data
  features_by_cells@meta.data$initial_clusters <- genes_by_cells@active.ident[match(rownames(features_by_cells@meta.data), names(genes_by_cells@active.ident))]
  plot_umap(features_by_cells, "umap_features_initial_clusters.jpg", color_by="initial_clusters")
  
}
