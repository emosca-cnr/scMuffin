#' landscent_sr
#' @param genes_by_cells real data
#' @param ppi The adjacency matrix of a user-given PPI network with rownames and colnames labeling genes (same gene identifier as in genes_by_cells)
#' @param reduceMethod indicates the method to do dimension reduction: "PCA" or "tSNE"
#' @param clusterMethod indicates the method to do clustering: "dbscan" or "PAM"
#' @return landscent_list list with the following elements: SR, DPT, potency_states, complete_output
#' @import LandSCENT destiny

landscent_sr <- function(genes_by_cells, ppi=NULL, reduceMethod = "PCA", clusterMethod = "dbscan", mc.cores=2){
	
	#integration matrix - ppi
	if(is.null(ppi)){
		cat("Using LandSCENT::net17Jan16.m\n")
		ppi <- scMuffin::ppi
	}
	
	temp <- sum(rownames(genes_by_cells) %in% rownames(ppi))
	if(temp<5){
		stop("less than 5 elements in common between geens-by-cells matrix and ppi. Check the identifiers!")
	}
	
	if(min(genes_by_cells)==0){
		cat("Detected 0 values...")
		half_min <- min(genes_by_cells[genes_by_cells>0]) / 2
		cat("\tsetting 0 values to:", half_min, "\n")
		genes_by_cells[genes_by_cells==0] <- half_min
	}
	Integration.l <- DoIntegPPI(exp.m = genes_by_cells, ppiA.m = ppi)
	
	#compute SR
	cat("Computing SR...\n")
	SR.o <- CompSRana(Integration.l, local = TRUE, mc.cores = mc.cores)
	
	#Infer the potency states
	cat("Inferring potency states...\n")
	InferPotency.o <- InferPotency(SR.o)
	PS <- InferPotency.o$potencyState
	
	#cat("Inferring landmark...\n")
	#InferLandmark.o <- InferLandmark(InferPotency.o, pheno.v = InferPotency.o$potencyState, reduceMethod = "PCA", clusterMethod = "dbscan")
	
	#`DoDiffusionMap` function
	cat("Diffusion map...\n")
	DoDiffusionMap.o <- DoDiffusionMap(InferPotency.o, mean_gap = 1, sd_gap = 1, root = c("cell", "state"), num_comp = 3, k = 30)
	
	#SR <- DoDiffusionMap.o$SR
	names(SR) <- colnames(DoDiffusionMap.o$data)
	
	dm <- DoDiffusionMap.o$DM
	root.idx <- IDoDiffusionMap.o$root
	
	cat("Diffusion pseudotimes...\n")
	dpt <- destiny::DPT(dm, tips = DoDiffusionMap.o$root)
	names(dpt) <- rownames(dm@eigenvectors)
	
	landscent_list <- list(SR= SR, DPT = dpt, potency_states = PS, complete_output = DoDiffusionMap.o)
	
	return(landscent_list)
	
}
