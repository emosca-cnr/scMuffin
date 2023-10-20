#' Create empirical null for gene set scoring
#' @param genes_by_cells genes-by-cells matrix
#' @param bins bins of genes, created by means of \code{\link{sc_data_bin}}
#' @param gene_set a vector of genes
#' @param k number of permutations
#' @importFrom stats na.omit
#' @description sc_create_null
#' @return A list of null gene sets.
#' @export

sc_create_null <- function(genes_by_cells=NULL, bins=NULL, gene_set=NULL, k=100){

  gene_set_idx <- na.omit(match(gene_set, rownames(genes_by_cells)))
  ans <- bins[gene_set_idx] ###bins that contain the gene set
  
  #remove the gene of the gene sets from the bins
  bins <- bins[-gene_set_idx]
  ans <- lapply(ans, function(ibin) which(bins == ibin)) #genes that belong the the bins (ibin) of the gene set

  kk <- min(unlist(lapply(ans, function(ibin) min(k, length(ibin))))) # the maximum actual k i constrained by ibin lenght, that is how many genes belong to such bin

  ans <- lapply(ans, function(ibin) sample(ibin, kk)) #sample a number of rows per bin
  ans <- lapply(ans, function(ibin) genes_by_cells[ibin, ]) #retrieve the corresponding data

  #old version: a unique df
  #ans <- do.call(rbind, ans)

  #new version: returns a list of null gene sets
  ans <- lapply(1:kk, function(x) do.call(rbind, lapply(ans, function(y) y[x, ])))

  return(ans)

}
