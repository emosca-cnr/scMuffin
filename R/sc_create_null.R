#' sc_create_null
#' @param seurat_data seurat data object
#' @param nbins number of bins
#' @param use.log use logarithm

sc_create_null <- function(seurat_data, bins, gene_set, k=100){

  ans <- bins[na.omit(match(gene_set, rownames(seurat_data)))]

  ans <- lapply(ans, function(ibin) which(bins == ibin)) #elements of bins that belong to ibin, i.e. rows of seurat data that belong to "ibin" set

  kk <- min(unlist(lapply(ans, function(ibin) min(k, length(ibin))))) # the maximum actual k i constrained by ibin lenght

  ans <- lapply(ans, function(ibin) sample(ibin, kk)) #sample a number of rows per bin
  ans <- lapply(ans, function(ibin) seurat_data[ibin, ]) #retrieve the corresponding data

  #old version: a unique df
  #ans <- do.call(rbind, ans)

  #new version: returns a list of null gene sets
  ans <- lapply(1:kk, function(x) do.call(rbind, lapply(ans, function(y) y[x, ])))

  return(ans)

}
