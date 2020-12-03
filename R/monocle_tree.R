#' monocle trajectory
#' 
#' @export
#' @import Seurat monocle DESeq2

monocle_tree <- function(genes_by_cells, root=NULL){
	
	cat("creating monocle cell data object...\n")
	
	fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name=rownames(genes_by_cells), row.names = rownames(genes_by_cells), stringsAsFactors = F))
	gbc_monocle <- monocle::newCellDataSet(genes_by_cells, featureData = fd, expressionFamily = VGAM::uninormal())
	
	
	cat("reducing dimensions...\n")
	gbc_monocle <- monocle::reduceDimension(gbc_monocle, norm_method = "none", pseudo_expr = 0, scaling = F, relative_expr = F)
	
	cat("ordering cells...\n")
	gbc_monocle <- monocle::orderCells(gbc_monocle)
	
	if(!is.null(root)){
		cat("re-ordering cells by root ", root, "...\n")
		root_state <- gbc_monocle@phenoData@data$State[rownames(gbc_monocle@phenoData@data) == root]
		gbc_monocle <- monocle::orderCells(gbc_monocle, root_state = root_state)
	}
	
	return(gbc_monocle)
	
}