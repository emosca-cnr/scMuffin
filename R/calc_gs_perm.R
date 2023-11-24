#' Calculate permutations
#' @param rll numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory
#' @param perm vector of permuted names
#' @param gs gene set
#' @description Calculate permutations

calc_gs_perm <- function(rll=NULL, perm=NULL, gs=NULL){

  out <- setNames(numeric(length(rll)), names(rll))
  
  for(i in 1:length(rll)){
    out[i] <- es(which(perm[[i]] %in% gs), array(rll[[i]], dimnames=list(perm[[i]])))[, 1]
  }
  
  return(out)

}
