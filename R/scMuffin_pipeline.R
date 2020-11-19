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
	
	res_signatures <- calculate_signatures(genes_by_cells, custom_signatures=custom_signatures, mc.cores=mc.cores)

	dir.create("signatures")
	save(res_signatures$signatures_by_cells, file="signatures/signatures_by_cells.RData", compress = "bzip2")
	
	#output
	heatmap_signatures(res_signatures$signatures_by_clusters, file="signatures/heatmap_signatures.jpg")
	write.table(res_signatures$signatures_by_clusters, file="signatures/signatures_by_clusters.txt", sep = "\t", row.names = T, col.names = NA)
	
	
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
	dir.create("monocle")

	mon_res <- monocle_tree(genes_by_cells)
	save(mon_res, file="monocle/mon_res.RData", compress="bzip2")
	
	jpeg("monocle/traj.jpg", width = 180, height = 180, res=300, units="mm")
	monocle::plot_cell_trajectory(mon_res)
	dev.off()
	
	cat("monocle cnv...\n")
	mon_res_cnv <- monocle_tree(cnv_res)
	save(mon_res_cnv, file="monocle/mon_res_cnv.RData", compress="bzip2")
	
	jpeg("monocle/traj_cnv.jpg", width = 180, height = 180, res=300, units="mm")
	monocle::plot_cell_trajectory(mon_res_cnv)
	dev.off()
	

	##################	merge everithing   ##################	
	dir.create("features")
	
	feature_list <- list(
		data.frame(id=colnames(res_signatures$signatures_by_cells), t(res_signatures$signatures_by_cells), stringsAsFactors = F),
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
	
	# features
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
	

	#ESTABLISH FEATURE TYPES
	feature_type <- unlist(lapply(cells_by_features_df, class))
	
	##QUANTITATIVE feature correlation
	feature_corr <- cor(as.matrix(cells_by_features_df[, feature_type == "numeric"]), method="spearman")
	heatmap_features_corr(feature_corr, file = "features/heatmap_corr.jpg")
	
	##QUALITATIVE FEATURES ACROSS GLOBAL EXPRESSION CLUSTERS
	cat("Enrichment of quantitative features across clusters")
	temp <- cluster_gsea(as.matrix(GetAssayData(features_by_cells)), genes_by_cells@active.ident[match(colnames(features_by_cells), names(genes_by_cells@active.ident))])
	write.table(temp$nes, file = paste0("expr_clusters/cluster_nes.txt"), row.names = T, col.names = NA, sep="\t")
	write.table(temp$fdrq, file = paste0("expr_clusters/cluster_nes_fdrq.txt"), row.names = T, col.names = NA, sep="\t")

	temp$fdrq[temp$nes<0] <- 1
	top_features <- lapply(split(temp$fdrq, rownames(temp$fdrq)), function(x) colnames(temp$fdrq)[rank(x) <= 3 & x < 0.05])
	dir.create("expr_clusters/boxplots_quantitative_features")
	boxplot_cluster(GetAssayData(features_by_cells, slot = "scale.data"), genes_by_cells@active.ident[match(colnames(features_by_cells), names(genes_by_cells@active.ident))], top_features = top_features, dir_out = "expr_clusters/boxplots_quantitative_features")

	
	##QUALITATIVE FEATURES ACROSS GLOBAL EXPRESSION CLUSTERS
	res_chisq <- cluster_chisq(cells_by_features_df[, feature_type != "numeric"], cell_clusters = genes_by_cells@active.ident[match(rownames(cells_by_features_df), names(genes_by_cells@active.ident))])
	top_features <- lapply(split(t(res_chisq$fdr), colnames(res_chisq$fdr)), function(x) rownames(res_chisq$fdr)[rank(x) <= 2 & x < 0.05])

	dotplot_cluster(cells_by_features_df[, fact_columns], genes_by_cells@active.ident[match(rownames(cells_by_features_df), names(genes_by_cells@active.ident))], dir_out = "expr_clusters/dotplots_qualitative_features", top_features = top_features, cont_tables = res_chisq$ct)

	write.table(res_chisq$fdr, file = "expr_clusters/cluster_chisq.txt", row.names = T, col.names = NA, sep="\t")
	for(i in 1:nrow(res_chisq$fdr)){
		for(j in 1:ncol(res_chisq$fdr)){
			n <- (i-1)*ncol(res_chisq$fdr)+j
			write.table(res_chisq$ct[[n]], file = paste0("expr_clusters/cluster_ct_", rownames(res_chisq$fdr)[i], "_", colnames(res_chisq$fdr)[j], "_hisq.txt"), row.names = T, col.names = NA, sep="\t")
		}
	}
	
	#####features-by-clustering
	clusterings <- split(t(cells_by_features_df[, feature_type != "numeric"]), colnames(cells_by_features_df)[feature_type != "numeric"])
	for(i in 1:length(clusterings)){
		clusterings[[i]] <- array(paste0(names(clusterings)[i], "_", clusterings[[i]]), dimnames = list(rownames(cells_by_features_df)))
	}
	clusterings$glob_expr <- genes_by_cells$seurat_clusters
	
	cluster_gsea_res_nes <- clusterings
	cluster_gsea_res_fdr <- clusterings
	res_signatures$signatures_by_cells[is.na(res_signatures$signatures_by_cells)] <- 0
	for(i in 1:length(clusterings)){
		cat(i)
		#cluster_gsea_res_nes[[i]] <- cluster_gsea(res_signatures$signatures_by_cells, clusterings[[i]])
		#cluster_gsea_res_nes[[i]] <- lapply(cluster_gsea_res_nes[[i]], function(x) data.frame(ID=colnames(x), t(x), stringsAsFactors = F))
		cluster_gsea_res_fdr[[i]] <- cluster_gsea_res_nes[[i]]$fdrq
		cluster_gsea_res_nes[[i]] <- cluster_gsea_res_nes[[i]]$nes
	}
	
	cluster_gsea_res_nes_merged <- Reduce(merge, cluster_gsea_res_nes)
	rownames(cluster_gsea_res_nes_merged) <- cluster_gsea_res_nes_merged$ID
	cluster_gsea_res_nes_merged$ID <- NULL
	
	cluster_gsea_res_fdr_merged <- Reduce(merge, cluster_gsea_res_fdr)
	rownames(cluster_gsea_res_fdr_merged) <- cluster_gsea_res_fdr_merged$ID
	cluster_gsea_res_fdr_merged$ID <- NULL
	
	save(cluster_gsea_res_nes_merged, cluster_gsea_res_fdr_merged, file="features/features_by_clusterings.RData", compress = "bzip2")
	
	
}
