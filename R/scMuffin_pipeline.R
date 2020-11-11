#' scMuffin: main function that runs the whole pipeline
#' @param genes_by_cells seurat_object_data, genes-by-cells input matrix
#' @param custom_signatures list of gene sets
#' @import parallel stats
#' @importFrom utils write.table
#' @export

scMuffin_pipeline <- function(genes_by_cells, custom_signatures=NULL, gene_sets=NULL, mc.cores=2, reference=NULL, features_for_reclustering=c("SIG")){
	
	
	#### global clusters ###
	dir.create("expr_clusters")
	pal <- rainbow(length(levels(genes_by_cells@active.ident)))
	plot_umap(genes_by_cells, "expr_clusters/umap_expr.jpg", pal = pal)
	
	################# GENE SETS FOR CLUSTERING #################
	if(!is.null(gene_sets)){
		
		dir.create("gene_sets_clusters")
		gene_sets_clusters <- calc_cluster_by_gs(genes_by_cells, gene_sets)
		
	}
	
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
	dir.create("features")
	
	feature_list <- list(
		data.frame(id=colnames(SC_signatures_by_cell_matrix), t(SC_signatures_by_cell_matrix), stringsAsFactors = F),
		data.frame(id=names(exp_rate_score), expr=exp_rate_score, stringsAsFactors = F),
		data.frame(id=rownames(output_landscent), output_landscent, stringsAsFactors = F),
		data.frame(id=names(heatmap_CNV_clusters), cnv=heatmap_CNV_clusters, stringsAsFactors = F),
		data.frame(id=colnames(mon_res), state=mon_res$State, pt=mon_res$Pseudotime, stringsAsFactors = F),
		data.frame(id=colnames(mon_res_cnv), state_cnv=mon_res_cnv$State, pt_cnv=mon_res_cnv$Pseudotime, stringsAsFactors = F)
		#data.frame(id=colnames(genes_by_cells), expr_clusters=genes_by_cells@active.ident, stringsAsFactors = F)
	)
	
	if(!is.null(gene_sets)){
		gene_sets_clusters <- lapply(gene_sets_clusters, function(x) data.frame(id=colnames(genes_by_cells), x, stringsAsFactors = F))
		for(i in 1:length(temp)){
			colnames(gene_sets_clusters) <- names(gene_sets_clusters)[i]
		}
		feature_list <- list(feature_list, gene_sets_clusters)
	}
	
	cells_by_features_df <- merge_matrix(feature_list)
	
	# RE-CLUSTERING on the basis of the features
	reclust_features <- sort(unique(unlist(lapply(features_for_reclustering, function(x) which(grepl(x, colnames(cells_by_features_df)))))))
	cat("reclustering on the basis of")
	print(colnames(cells_by_features_df)[reclust_features])
	
	features_by_cells <- cells_by_features_df[, reclust_features]
	features_by_cells[is.na(features_by_cells)] <- 0
	features_by_cells <- re_clustering(t(as.matrix(features_by_cells)), n_comp = 2)
	save(features_by_cells, file="features/features_by_cells.RData", compress = "bzip2")
	
	#adding the reclustering info
	cells_by_features_df$feature_clusters <- features_by_cells@active.ident[match(rownames(cells_by_features_df), colnames(features_by_cells))]
	
	save(cells_by_features_df, file="features/cells_by_features_df.RData", compress = "bzip2")
	
	####PLOT UMAP BY EXPRESSION COLORED BY EVERY FEATURE
	plot_umap_expr_features(seurat_object, cells_by_features_df, dir="expr_clusters/")
	
	####PLOT UMAP BY FEATURES COLORED BY EVERY FEATURE
	plot_umap_expr_features(features_by_cells, cells_by_features_df, dir="features/")

	
	#ESTABLISH FEATURE TYPES
	feature_type <- unlist(lapply(cells_by_features_df, class))
	
	##QUANTITATIVE feature correlation
	feature_corr <- cor(as.matrix(cells_by_features_df[, feature_type == "numeric"]), method="spearman")
	heatmap_features_corr(feature_corr, file = "features/heatmap_corr.jpg")
	
	#trend of quantitative features in expression clusteres
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), genes_by_cells@active.ident, dir_out = "clusters/boxplots_quantitative_features")
	
	##QUALITATIVE
	feature_overlap <- table(cells_by_features_df[, feature_type != "numeric"])
	print(feature_overlap)
	feature_overlap <- as.data.frame(feature_overlap)
	feature_overlap <- feature_overlap[order(-feature_overlap$Freq), ]
	feature_overlap <- feature_overlap[feature_overlap$Freq > 0, ]
	write.table(feature_overlap, file="features/qualitative_features_overlap.txt", row.names = T, col.names = NA, sep="\t")
	
}
