#' Computation of Diffusion Map and Diffusion Pseudo Time
#' 
#' @param scMuffinList scMuffinList object
#' @param root_cell root cell or "random" to pick up a random root via the destiny::random_root function
#' @param n_pcs number of PCs to consider, 50 by default
#' @param ... parameters passed to destiny::DiffusionMap
#' @importFrom destiny DiffusionMap random_root DPT
#' @description The functions computes a Diffusion Map and Diffusion Pseudo Time using package destiny.
#' @return scMuffinList with element "diffusion_map_pseudo_t", a list with summary and full. summary is a data.frame with: the first two eigenvectors; branch, branch label for each cell; tips, whether the cell is a tip of the branching. full contains dpt object, generated by [destiny::DPT()]. See destiny for further information.
#' @export

diff_map <- function(scMuffinList = NULL, root_cell="random", n_pcs=50, ...){
  
  cat("Calculating the diffusion map...\n")
  res <- DiffusionMap(t(as.matrix(scMuffinList$normalized)), n_pcs=n_pcs, ...)
  
  #Il DPT ?? metrica che dipende dalla cellula scelta, identical(dpt[root.idx], dpt$dpt) dpt[["dpt"]]
  dpt <- NA
  
  if(root_cell == "random"){
    cat("calculating a random root...\n")
    root_cell <-random_root(res) #NUll
  }
  cat("Diffusion pseudotimes...\n")
  dpt <- DPT(res, tips = root_cell)

  scMuffinList$diffusion_map_pseudo_t <- list(
    summary = data.frame(res@eigenvectors[, 1:2], dpt=dpt$dpt, branch=dpt@branch[, 1], tips=dpt@tips[, 1], stringsAsFactors = F, row.names = rownames(res@eigenvectors)),
    full = dpt
  )
  
  return(scMuffinList)
  
}
