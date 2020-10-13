#' landscent
#' @param genes_by_cells_matrix real data
#' @param ppi The adjacency matrix of a user-given PPI network with rownames and colnames labeling genes (same gene identifier as in genes_by_cells_matrix)
#' @param reduceMethod indicates the method to do dimension reduction: "PCA" or "tSNE"
#' @param clusterMethod indicates the method to do clustering: "dbscan" or "PAM"
#' @return landscent_list list with the following elements: SR, DPT, potency_states, complete_output
#' @import LandSCENT
#' @importFrom destiny DPT

landscent <- function(genes_by_cells_matrix_matrix, ppi=NULL, reduceMethod = "PCA", clusterMethod = "dbscan", mc.cores=2){
	
	
	#integration matrix - ppi
	if(is.null(ppi)){
		cat("Using LandSCENT::net17Jan16.m\n")
		ppi <- scMuffin::ppi
	}
	
	temp <- sum(rownames(genes_by_cells_matrix) %in% rownames(ppi))
	if(temp<5){
		stop("less than 5 elements in common between geens-by-cells matrix and ppi. Check the identifiers!")
	}
	
	if(min(genes_by_cells_matrix)==0){
		cat("Detected 0 values...")
		half_min <- min(genes_by_cells_matrix[genes_by_cells_matrix>0]) / 2
		cat("\tsetting 0 values to:", half_min, "\n")
		genes_by_cells_matrix[genes_by_cells_matrix==0] <- half_min
	}
	Integration.l <- LandSCENT::DoIntegPPI(exp.m = genes_by_cells_matrix, ppiA.m = ppi)
	
	#compute Signalling entRopy
	cat("Computing SR...\n")
	Integration.l <- LandSCENT::CompSRana(Integration.l, local = TRUE, mc.cores = mc.cores)
	
	#Infer the potency states
	cat("Inferring potency states...\n")
	Integration.l <- LandSCENT::InferPotency(Integration.l)
	
	#cat("Inferring landmark...\n")
	#InferLandmark.o <- InferLandmark(InferPotency.o, pheno.v = InferPotency.o$potencyState, reduceMethod = "PCA", clusterMethod = "dbscan")
	
	#extract DiffusionMap and root index objects to calculate DPT
	cat("Diffusion map...\n")
	Integration.l <- LandSCENT::DoDiffusionMap(Integration.l, mean_gap = 1, sd_gap = 1, root = c("cell", "state"), num_comp = 3, k = 30)
	
	#Il DPT è metrica che dipende dalla cellula scelta, identical(dpt[root.idx], dpt$dpt) dpt[["dpt"]]
	cat("Diffusion pseudotimes...\n")
	dpt <- destiny::DPT(Integration.l$DM, tips = Integration.l$root)
	
	#Create DF: cell name, dpt, SR and PS
	score <- data.frame(dpt=dpt$dpt, SR=Integration.l$SR, PS=Integration.l$potencyState, stringsAsFactors = F, row.names = rownames(Integration.l$DM@eigenvectors))
	
	return(score)
	
}
