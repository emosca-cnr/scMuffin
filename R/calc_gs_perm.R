#' Calculate permutations
#' @param rll numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory
#' @param perm vector of permuted names
#' @param gs gene set
#' @description Calculate permutations

calc_gs_perm <- function(rll, perm, gs){

  #out <- es(which(perm %in% gs), array(rl, dimnames=list(perm)))
  out <- unlist(lapply(rll, function(x) es(which(perm %in% gs), array(x, dimnames=list(perm)))))

  return(out)

}
