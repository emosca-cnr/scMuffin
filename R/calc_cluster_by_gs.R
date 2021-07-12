#' calc_cluster_by_gs
#'
#'
#' @export


calc_cluster_by_gs <- function(genes_by_cells, gene_sets, out_dir="./"){
	
	if(!dir.exists(out_dir)){
		dir.create(out_dir)
	}
	
	gene_sets_clusters <- vector("list", length(gene_sets))
	for(i in 1:length(gene_sets)){
		
		res <- cluster_by_gs(genes_by_cells, gs = gene_sets[[i]])
		gene_sets_clusters[[i]] <- res$seurat_clusters
		plot_umap(res, file = paste0(out_dir, "/umap_", names(gene_sets)[i], ".jpg"), group.by="seurat_clusters")
		
	}
	
	names(gene_sets_clusters) <- names(gene_sets)
	
	gene_sets_clusters <- lapply(gene_sets_clusters, function(x) data.frame(id=colnames(genes_by_cells), x, stringsAsFactors = F))
	for(i in 1:length(gene_sets_clusters)){
		colnames(gene_sets_clusters[[i]]) <- c("id", names(gene_sets_clusters)[i])
	}
	
	gene_sets_clusters <- Reduce(merge, gene_sets_clusters)
	rownames(gene_sets_clusters) <- gene_sets_clusters$id
	gene_sets_clusters$id <- NULL
	
	return(gene_sets_clusters)
}
