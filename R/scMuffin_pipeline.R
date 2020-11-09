#' scMuffin: main function that runs the whole pipeline
#' @param genes_by_cells seurat_object_data, genes-by-cells input matrix
#' @param custom_signatures list of gene sets
#' @import parallel stats
#' @importFrom utils write.table
#' @export

scMuffin_pipeline <- function(genes_by_cells, custom_signatures=NULL, mc.cores=2, reference=NULL){
	
	
	dir.create("clusters")
	
	#### Visual comparison between initial and final clusters ###
	pal <- rainbow(length(levels(genes_by_cells@active.ident)))
	plot_umap(genes_by_cells, "clusters/umap_by_gene_expression.jpg", pal = pal)
	
	
	##################	SIGNATURES #################
	cat("Calcuting gene signature scores...\n")
	
	SC_signatures_by_cluster_matrix <- calculate_signatures(genes_by_cells, custom_signatures=custom_signatures, mc.cores=mc.cores)
	
	#output
	heatmap_signatures(SC_signatures_by_cluster_matrix, file="signatures/heatmap_signatures.jpg")
	write.table(SC_signatures_by_cluster_matrix, file="signatures/signatures_by_clusters.txt", sep = "\t", row.names = T, col.names = NA)
	
	

	##################	EXPRESSION RATE   ##################	
	exp_rate_score <- exp_rate(as.matrix(genes_by_cells@assays$RNA@counts))
	dir.create("expr_score")
	save(exp_rate_score, file="expr_score/exp_rate_score.RData", compress = "bzip2")
	

	##################	CNV @Valentina   ##################	
	
	cnv_res <- calculate_CNV(as.data.frame(GetAssayData(genes_by_cells)), mc.cores = mc.cores, reference = reference)
	
	ngenes_chrom <- unlist(lapply(cnv_res, nrow)) # number of genes per chromosome
	cnv_res <- preprocess_for_heatmap_CNV(cnv_res)
	
	dir.create("cnv")
	heatmap_CNV_clusters <- heatmap_CNV(cnv_res, ngenes_chrom, file = "cnv/heatmap_CNV.jpg", reference = "reference")
	save(cnv_res, heatmap_CNV_clusters, file="cnv/cnv_res.RData", compress = "bzip2")
	
	genes_by_cells@meta.data$cnv <- heatmap_CNV_clusters[match(rownames(genes_by_cells@meta.data), names(heatmap_CNV_clusters))]
	pal <- rainbow(length(levels(genes_by_cells@meta.data$cnv)))
	plot_umap(genes_by_cells, "clusters/umap_genes_cnv_clusters.jpg", color_by="cnv", pal = pal)
	
	##################	Signalling entropy rate (SR) @Noemi 	  ##################	
	##################	Potency states (LandScent): labels @Noemi 	  ##################	
	##################	Diffusion pseudotime (DPT) (Destiny): @Noemi   ##################	
	
	cat("Calcuting landscent-related scores...\n")
	output_landscent <- landscent(as.matrix(genes_by_cells@assays$RNA@data), mc.cores=mc.cores)
	
	dir.create("landscent")
	save(output_landscent, file="landscent/output_landscent.RData", compress = "bzip2")
	
	################## Monocle ###########
	cat("monocle...\n")
	mon_res <- monocle_tree(genes_by_cells)
	dir.create("monocle")
	save(mon_res, file="monocle/mon_res.RData", compress="bzip2")
	
	jpeg("monocle/traj.jpg", width = 180, height = 180, res=300, units="mm")
	monocle::plot_cell_trajectory(mon_res)
	dev.off()
	
	cat("monocle cnv...\n")
	mon_res_cnv <- monocle_tree(cnv_res)
	dir.create("monocle")
	save(mon_res_cnv, file="monocle/mon_res_cnv.RData", compress="bzip2")
	
	jpeg("monocle/traj_cnv.jpg", width = 180, height = 180, res=300, units="mm")
	monocle::plot_cell_trajectory(mon_res_cnv)
	dev.off()
	
	
	##################	merge everithing   ##################	
	feature_list <- list(
		data.frame(id=colnames(SC_signatures_by_cell_matrix), t(SC_signatures_by_cell_matrix), stringsAsFactors = F),
		data.frame(id=names(exp_rate_score), expr=exp_rate_score, stringsAsFactors = F),
		data.frame(id=rownames(output_landscent), output_landscent, stringsAsFactors = F),
		data.frame(id=names(heatmap_CNV_clusters), cnv=heatmap_CNV_clusters, stringsAsFactors = F),
		data.frame(id=colnames(mon_res), state=mon_res$State, pt=mon_res$Pseudotime, stringsAsFactors = F),
		data.frame(id=colnames(mon_res_cnv), state_cnv=mon_res_cnv$State, pt_cnv=mon_res_cnv$Pseudotime, stringsAsFactors = F)
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
	
	heatmap_cluster_enrichment(clust_enr_res$nes, clust_enr_res$fdrq, file = "cluster_comparison/heatmap_cluster_enrichment_features.jpg")
	heatmap_cluster_enrichment(clust_enr_res_global_expr$nes, clust_enr_res_global_expr$fdrq, file = "cluster_comparison/heatmap_cluster_enrichment_glob_expr.jpg", cex.axis = 0.6)
	
	
	pal <- rainbow(length(levels(features_by_cells@active.ident)))
	plot_umap(features_by_cells, "clusters/umap_by_features.jpg", pal = pal)
	
	
	### COLOR BY FEATURE
	plot_umap_colored_features(genes_by_cells, features_by_cells, dir="clusters") 
	
	#Visualization of CNV on gene expression clusters
	
	features_by_cells@meta.data$cnv <- heatmap_CNV_clusters[match(rownames(features_by_cells@meta.data), names(heatmap_CNV_clusters))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$cnv)))
	plot_umap(features_by_cells, "clusters/umap_features_cnv_clusters.jpg", color_by="cnv", pal = pal)
	
	#add initial clusters information in features_by_cells meta.data
	features_by_cells@meta.data$initial_clusters <- genes_by_cells@active.ident[match(rownames(features_by_cells@meta.data), names(genes_by_cells@active.ident))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$initial_clusters)))
	plot_umap(features_by_cells, "clusters/umap_features_initial_clusters.jpg", color_by="initial_clusters", pal = pal)
	
	#potency state plot
	genes_by_cells@meta.data$ps <- factor(output_landscent$PS[match(rownames(genes_by_cells@meta.data), rownames(output_landscent))])
	pal <- rainbow(length(levels(genes_by_cells@meta.data$ps)))
	plot_umap(genes_by_cells, "clusters/umap_genes_ps_clusters.jpg", color_by="ps", pal = pal)
	
	features_by_cells@meta.data$ps <- factor(output_landscent$PS[match(rownames(features_by_cells@meta.data), rownames(output_landscent))])
	pal <- rainbow(length(levels(features_by_cells@meta.data$ps)))
	plot_umap(features_by_cells, "clusters/umap_features_ps_clusters.jpg", color_by="ps", pal = pal)
	
	#BOXPLOT feature in clusters
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), genes_by_cells@active.ident, dir_out = "clusters/boxplot_cluster_by_genes")
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), features_by_cells@active.ident, dir_out = "clusters/boxplot_cluster_by_features")
	
}
