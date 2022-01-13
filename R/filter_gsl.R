#' Filter a gene set list
#' @param gsl gene set list (named list)
#' @param universe set of all possible values for items of gsl
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @export
#' 
filter_gsl <- function(gsl, universe, min_size=5, max_size=500){
  
  ans <- lapply(gsl, function(x) x[x %in% universe])
  idx_keep <- unlist(lapply(ans, function(x) length(x) >= min_size & length(x) <= max_size))
  
  ans <- ans[idx_keep]
  
  return(ans)
  
}