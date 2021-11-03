#' Computation of Diffusion Map
#' 
#' The function apply Landscent package to a gene by cell matric to compute signal entropy, potency states and diffusion pseudo time
#' @param Seurat_object Seurat object with normalized data
#' @param output_ndc The numnber of diffusion components that will be included in the output
#' @param root_cell root cell or "random" to pick up a random root via the destiny::random_root function
#' @param ... parameters for destiny::DiffusionMap
#' @importFrom destiny DiffusionMap random_root DPT
#' @export

diff_map <- function(Seurat_object, output_ndc=1:3, root_cell=NULL, ...){
	
  cat("Calculating the diffusion map...\n")
  res <- destiny::DiffusionMap(t(as.matrix(GetAssayData(Seurat_object, slot = "data"))), ...)
	
  #Il DPT ?? metrica che dipende dalla cellula scelta, identical(dpt[root.idx], dpt$dpt) dpt[["dpt"]]
  dpt_raw <- dpt <- NA
  if(!is.null(root_cell)){
    if(root_cell == "random"){
      cat("calculating a random root...\n")
      root_cell <- destiny::random_root(res) #NUll
    }
	  cat("Diffusion pseudotimes...\n")
	  dpt <- destiny::DPT(res, tips = root_cell)
	  dpt_raw <- dpt$dpt
	  dpt <- -dpt$dpt + max(dpt$dpt)
	}
	
	#Create DF: cell name, dpt, SR and PS
	score <- data.frame(res@eigenvectors[, output_ndc], dpt_raw=dpt_raw, dpt=dpt, stringsAsFactors = F, row.names = rownames(res@eigenvectors))
	
	return(score)
	
}
