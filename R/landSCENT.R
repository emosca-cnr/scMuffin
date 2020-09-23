#' LandSCENT
#' @param genes_by_cells real data
#' @param ppi The adjacency matrix of a user-given PPI network with rownames and colnames labeling genes (same gene identifier as in genes_by_cells)
#' @param reduceMethod indicates the method to do dimension reduction: "PCA" or "tSNE"
#' @param clusterMethod indicates the method to do clustering: "dbscan" or "PAM"
#' @return landscent_list list with the following elements: SR, DPT, potency_states, complete_output

landSCENT <- function(genes_by_cells, ppi, reduceMethod = "PCA", clusterMethod = "dbscan", ){
  
#integration matrix - ppi
Integration.l <- DoIntegPPI(exp.m = genes_by_cells, ppiA.m = ppi, log_trans = T)

#compute SR
SR.o <- CompSRana(Integration.l, local = TRUE, mc.cores = 4)

#Infer the potency states
InferPotency.o <- InferPotency(SR.o)
PS <- InferPotency.o$potencyState

InferLandmark.o <- InferLandmark(InferPotency.o, pheno.v = InferPotency.o$potencyState,
                                 reduceMethod = "PCA", clusterMethod = "dbscan")

#`DoDiffusionMap` function
DoDiffusionMap.o <- DoDiffusionMap(InferPotency.o,
                                   mean_gap = 1, sd_gap = 1,
                                   root = c("cell", "state"),
                                   num_comp = 3,
                                   k = 30)

SR <- DoDiffusionMap.o$SR
dm <- DoDiffusionMap.o$DM
root.idx <- IDoDiffusionMap.o$root
dpt <- destiny::DPT(dm, tips = DoDiffusionMap.o$root)

landscent_list <- list(SR= SR,
                       DPT = dpt,
                       potency_states = PS,
                       complete_output = DoDiffusionMap.o)
return(landscent_list)
}
