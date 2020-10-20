#' scMuffin: main function that runs the whole pipeline
#' @param genes_by_cells seurat_object_data, genes-by-cells input matrix
#' @param custom_signatures list of gene sets
#' @import parallel stats
#' @importFrom utils write.table
#' @export

scMuffin_pipeline <- function(genes_by_cells, custom_signatures=NULL, analyses = c("signatures", "landscent", "cnv", "expr_rate"), mc.cores=2){
	
	
	##################	SIGNATURES #################
	SC_signatures_by_cell_matrix <- NULL
	if("signatures" %in% analyses){
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
		
		save(SC_signatures_by_cell_matrix, file="SC_signatures_by_cell_matrix.RData", compress = "bzip2")
		
		#gene-set score per cluster list
		res_signatures_clusters <- lapply(res_signatures, function(i_marker_res) gene_set_score_in_clusters(i_marker_res$score_table, genes_by_cells@active.ident, ncells_min = 5))
		
		#signatures-by-clusters matrix
		SC_signatures_by_cluster_matrix <- do.call(rbind, lapply(res_signatures_clusters, function(x) array(x$score[order(x$cluster)], dimnames = list(c(x$cluster[order(x$cluster)])))))
		
		#output
		heatmap_signatures(SC_signatures_by_cluster_matrix)
		write.table(SC_signatures_by_cluster_matrix, file="signatures_by_clusters.txt", sep = "\t", row.names = T, col.names = NA)
		
	}
	##################	Signalling entropy rate (SR) @Noemi 	  ##################	
	##################	Potency states (LandScent): labels @Noemi 	  ##################	
	##################	Diffusion pseudotime (DPT) (Destiny): @Noemi   ##################	
	output_landscent <- NULL
	if("landscent" %in% analyses){
		
		cat("Calcuting landscent-related scores...\n")
		output_landscent <- landscent(as.matrix(genes_by_cells@assays$RNA@data), mc.cores=mc.cores)
		save(output_landscent, file="output_landscent.RData", compress = "bzip2")
		
	}
	
	##################	EXPRESSION RATE   ##################	
	exp_rate_score <- NULL
	if("expr_rate" %in% analyses){
		exp_rate_score <- exp_rate(as.matrix(genes_by_cells@assays$RNA@counts))
		save(exp_rate_score, file="exp_rate_score.RData", compress = "bzip2")
	}
	##################	Cell cycle state   ##################	
	
	
	##################	CNV @Valentina   ##################	
	heatmap_CNV_clusters <- NULL
	if("cnv" %in% analyses){
		cnv_res <- calculate_CNV(as.matrix(genes_by_cells@assays$RNA@data), mc.cores = mc.cores)
		
		ngenes_chrom <- unlist(lapply(cnv_res, nrow)) # number of genes per chromosome
		cnv_res <- preprocess_for_heatmap_CNV(cnv_res)
		heatmap_CNV_clusters <- heatmap_CNV(cnv_res, ngenes_chrom)
		save(cnv_res, heatmap_CNV_clusters, file="cnv_res.RData", compress = "bzip2")
	}
	
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
	heatmap_features_corr(feature_corr)
	
	##################	re-clustering   ##################	
	features_by_cells <- re_clustering(t(as.matrix(cells_by_features_df[, feature_type == "numeric"])))
	save(features_by_cells, file="features_by_cells.RData", compress = "bzip2")
	
	dendrogram_genes <- as.dendrogram(BuildClusterTree(genes_by_cells, reorder = T, features = rownames(genes_by_cells))@tools$BuildClusterTree)
	dendrogram_features <- as.dendrogram(BuildClusterTree(features_by_cells, reorder = T, features = rownames(features_by_cells), slot = "scale.data")@tools$BuildClusterTree)
	
	clust_comp_matrix <- cluster_comparison(genes_by_cells@active.ident, features_by_cells@active.ident, dendrogram_genes, dendrogram_features)
	
	plot_cluster_comparison(clust_comp_matrix, dend_row = dendrogram_genes, dend_col = dendrogram_features)
	
	#### FINAL OUTPUTS ###
	
	#### Visual comparison between initial and final clusters ###
	pal <- rainbow(length(levels(genes_by_cells@active.ident)))
	plot_umap(genes_by_cells, "umap_by_gene_expression.jpg", pal = pal)
	
	pal <- rainbow(length(levels(features_by_cells@active.ident)))
	plot_umap(features_by_cells, "umap_by_features.jpg", pal = pal)
	
	
	### COLOR BY FEATURE
	plot_umap_colored_features(genes_by_cells, features_by_cells)

	#Visualization of CNV on gene expression clusters
	genes_by_cells@meta.data$cnv <- heatmap_CNV_clusters[match(rownames(genes_by_cells@meta.data), names(heatmap_CNV_clusters))]
	pal <- rainbow(length(levels(genes_by_cells@meta.data$cnv)))
	plot_umap(genes_by_cells, "umap_genes_cnv_clusters.jpg", color_by="cnv", pal = pal)
	
	features_by_cells@meta.data$cnv <- heatmap_CNV_clusters[match(rownames(features_by_cells@meta.data), names(heatmap_CNV_clusters))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$cnv)))
	plot_umap(features_by_cells, "umap_features_cnv_clusters.jpg", color_by="cnv", pal = pal)
	
	#add initial clusters information in features_by_cells meta.data
	features_by_cells@meta.data$initial_clusters <- genes_by_cells@active.ident[match(rownames(features_by_cells@meta.data), names(genes_by_cells@active.ident))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$initial_clusters)))
	plot_umap(features_by_cells, "umap_features_initial_clusters.jpg", color_by="initial_clusters", pal = pal)
	
	#potency state plot
	genes_by_cells@meta.data$ps <- output_landscent$PS[match(rownames(genes_by_cells@meta.data), rownames(output_landscent))]
	pal <- rainbow(length(levels(genes_by_cells@meta.data$ps)))
	plot_umap(genes_by_cells, "umap_genes_ps_clusters.jpg", color_by="ps", pal = pal)
	
	features_by_cells@meta.data$ps <- output_landscent$PS[match(rownames(features_by_cells@meta.data), rownames(output_landscent))]
	pal <- rainbow(length(levels(features_by_cells@meta.data$ps)))
	plot_umap(features_by_cells, "umap_features_ps_clusters.jpg", color_by="ps", pal = pal)
	
	#BOXPLOT feature in clusters
	boxplot_cluster(GetAssayData(genes_by_cells, slot = "scale.data"), genes_by_cells@active.ident, dir_out = "boxplot_cluster_by_genes")
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), features_by_cells@active.ident, dir_out = "boxplot_cluster_by_features")
	
	
}
