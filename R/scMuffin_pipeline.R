#' scMuffin: main function that runs the whole pipeline
#' @param genes_by_cells seurat_object_data, genes-by-cells input matrix
#' @param custom_signatures list of gene sets
#' @import parallel stats
#' @importFrom utils write.table
#' @export

scMuffin_pipeline <- function(genes_by_cells, custom_signatures=NULL, mc.cores=2, reference=NULL){
	
	
	##################	SIGNATURES #################
	cat("Calcuting gene signature scores...\n")
	
	if(!is.null(custom_signatures)){
		signatures <- c(signatures, custom_signatures)
	}
	cat("# of signatures: ", length(signatures), "\n")
	
	#dataset bins
	data_bins <- sc_data_bin(as.matrix(GetAssayData(genes_by_cells)), nbins = 25, use.log = TRUE)
	
	#signatures-by-cell matrix
	res_signatures <- parallel::mclapply(signatures, function(i_marker_set) gene_set_score(i_marker_set, genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=100, nmark_min = 5, ncells_min = 5), mc.cores = mc.cores)
	
	SC_signatures_by_cell_matrix <- t(do.call(cbind, lapply(res_signatures, function(x) array(x$score_table$avg_delta_score, dimnames = list(rownames(x$score_table))))))
	
	dir.create("signatures")
	save(SC_signatures_by_cell_matrix, file="signatures/SC_signatures_by_cell_matrix.RData", compress = "bzip2")
	
	#gene-set score per cluster list
	res_signatures_clusters <- lapply(res_signatures, function(i_marker_res) gene_set_score_in_clusters(i_marker_res$score_table, genes_by_cells@active.ident, ncells_min = 5))
	
	#signatures-by-clusters matrix
	SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
	
	#output
	heatmap_signatures(SC_signatures_by_cluster_matrix, file="signatures/heatmap_signatures.jpg")
	write.table(SC_signatures_by_cluster_matrix, file="signatures/signatures_by_clusters.txt", sep = "\t", row.names = T, col.names = NA)
	
	
	##################	Signalling entropy rate (SR) @Noemi 	  ##################	
	##################	Potency states (LandScent): labels @Noemi 	  ##################	
	##################	Diffusion pseudotime (DPT) (Destiny): @Noemi   ##################	
	
	cat("Calcuting landscent-related scores...\n")
	output_landscent <- landscent(as.matrix(genes_by_cells@assays$RNA@data), mc.cores=mc.cores)
	
	dir.create("landscent")
	save(output_landscent, file="landscent/output_landscent.RData", compress = "bzip2")
	
	
	
	##################	EXPRESSION RATE   ##################	
	exp_rate_score <- exp_rate(as.matrix(genes_by_cells@assays$RNA@counts))
	dir.create("expr_score")
	save(exp_rate_score, file="expr_score/exp_rate_score.RData", compress = "bzip2")
	
	##################	Cell cycle state   ##################	
	
	
	##################	CNV @Valentina   ##################	
	
	cnv_res <- calculate_CNV(as.matrix(genes_by_cells@assays$RNA@data), mc.cores = mc.cores, reference = reference)
	
	ngenes_chrom <- unlist(lapply(cnv_res, nrow)) # number of genes per chromosome
	cnv_res <- preprocess_for_heatmap_CNV(cnv_res)
	
	dir.create("cnv")
	heatmap_CNV_clusters <- heatmap_CNV(cnv_res, ngenes_chrom, file = "cnv/heatmap_CNV.jpg", reference = "reference")
	save(cnv_res, heatmap_CNV_clusters, file="cnv/cnv_res.RData", compress = "bzip2")
	
	
	##################	merge everithing   ##################	
	feature_list <- list(
		data.frame(id=colnames(SC_signatures_by_cell_matrix), t(SC_signatures_by_cell_matrix), stringsAsFactors = F),
		data.frame(id=names(exp_rate_score), expr=exp_rate_score, stringsAsFactors = F),
		data.frame(id=rownames(output_landscent), output_landscent, stringsAsFactors = F),
		data.frame(id=names(heatmap_CNV_clusters), cnv=heatmap_CNV_clusters, stringsAsFactors = F)
	)
	
	cells_by_features_df <- merge_matrix(feature_list)
	
	feature_type <- unlist(lapply(cells_by_features_df, class))
	
	feature_corr <- cor(as.matrix(cells_by_features_df[, feature_type == "numeric"]), method="spearman")
	
	dir.create("features")
	heatmap_features_corr(feature_corr, file = "features/heatmap_corr.jpg")
	
	##################	re-clustering   ##################	
	features_by_cells <- re_clustering(t(as.matrix(cells_by_features_df[, feature_type == "numeric"])))
	save(features_by_cells, file="features/features_by_cells.RData", compress = "bzip2")
	
	dendrogram_genes <- as.dendrogram(BuildClusterTree(genes_by_cells, reorder = T, features = rownames(genes_by_cells))@tools$BuildClusterTree)
	dendrogram_features <- as.dendrogram(BuildClusterTree(features_by_cells, reorder = T, features = rownames(features_by_cells), slot = "scale.data")@tools$BuildClusterTree)
	
	clust_comp_matrix <- cluster_comparison(genes_by_cells@active.ident, features_by_cells@active.ident, dendrogram_genes, dendrogram_features)
	dir.create("cluster_comparison")
	plot_cluster_comparison(clust_comp_matrix$freq, dend_row = dendrogram_genes, dend_col = dendrogram_features, file = "cluster_comparison/heatmpa_cluster_comparison.jpg")
	
	clust_enr_res <- cluster_gsea(as.matrix(GetAssayData(features_by_cells)), features_by_cells@active.ident)
	write.table(clust_enr_res$nes, file="cluster_comparison/clust_enr_res_nes.txt", row.names = T, col.names = NA, sep="\t")
	write.table(clust_enr_res$fdrq, file="cluster_comparison/clust_enr_res_fdrq.txt", row.names = T, col.names = NA, sep="\t")
	
	clust_enr_res_global_expr <- cluster_gsea(as.matrix(GetAssayData(features_by_cells)), genes_by_cells@active.ident)
	write.table(clust_enr_res_global_expr$nes, file="cluster_comparison/clust_enr_res_glob_expr_nes.txt", row.names = T, col.names = NA, sep="\t")
	write.table(clust_enr_res_global_expr$fdrq, file="cluster_comparison/clust_enr_res_glob_expr_fdrq.txt", row.names = T, col.names = NA, sep="\t")
	
	#clust_chi_res <- cluster_chisq(as.matrix(GetAssayData(features_by_cells)), features_by_cells@active.ident)
	
	heatmap_cluster_enrichment(clust_enr_res$nes, clust_enr_res$fdrq, file = "cluster_comparison/heatmap_cluster_enrichment_features.jpg")
	heatmap_cluster_enrichment(clust_enr_res_global_expr$nes, clust_enr_res_global_expr$fdrq, file = "cluster_comparison/heatmap_cluster_enrichment_glob_expr.jpg", cex.axis = 0.6)
	
	
	#### FINAL OUTPUTS ###
	
	dir.create("clusters")
	
	#### Visual comparison between initial and final clusters ###
	pal <- rainbow(length(levels(genes_by_cells@active.ident)))
	plot_umap(genes_by_cells, "clusters/umap_by_gene_expression.jpg", pal = pal)
	
	pal <- rainbow(length(levels(features_by_cells@active.ident)))
	plot_umap(features_by_cells, "clusters/umap_by_features.jpg", pal = pal)
	
	
	### COLOR BY FEATURE
	plot_umap_colored_features(genes_by_cells, features_by_cells, dir="clusters") 
	
	#Visualization of CNV on gene expression clusters
	genes_by_cells@meta.data$cnv <- heatmap_CNV_clusters[match(rownames(genes_by_cells@meta.data), names(heatmap_CNV_clusters))]
	pal <- rainbow(length(levels(genes_by_cells@meta.data$cnv)))
	plot_umap(genes_by_cells, "clusters/umap_genes_cnv_clusters.jpg", color_by="cnv", pal = pal)
	
	features_by_cells@meta.data$cnv <- heatmap_CNV_clusters[match(rownames(features_by_cells@meta.data), names(heatmap_CNV_clusters))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$cnv)))
	plot_umap(features_by_cells, "clusters/umap_features_cnv_clusters.jpg", color_by="cnv", pal = pal)
	
	#add initial clusters information in features_by_cells meta.data
	features_by_cells@meta.data$initial_clusters <- genes_by_cells@active.ident[match(rownames(features_by_cells@meta.data), names(genes_by_cells@active.ident))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$initial_clusters)))
	plot_umap(features_by_cells, "clusters/umap_features_initial_clusters.jpg", color_by="initial_clusters", pal = pal)
	
	#potency state plot
	genes_by_cells@meta.data$ps <- output_landscent$PS[match(rownames(genes_by_cells@meta.data), rownames(output_landscent))]
	pal <- rainbow(length(levels(genes_by_cells@meta.data$ps)))
	plot_umap(genes_by_cells, "clusters/umap_genes_ps_clusters.jpg", color_by="ps", pal = pal)
	
	features_by_cells@meta.data$ps <- output_landscent$PS[match(rownames(features_by_cells@meta.data), rownames(output_landscent))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$ps)))
	plot_umap(features_by_cells, "clusters/umap_features_ps_clusters.jpg", color_by="ps", pal = pal)
	
	#BOXPLOT feature in clusters
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), genes_by_cells@active.ident, dir_out = "clusters/boxplot_cluster_by_genes")
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), features_by_cells@active.ident, dir_out = "clusters/boxplot_cluster_by_features")
	
	
}
