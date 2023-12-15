#' Calculate permutations
#' @param rll numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory
#' @param perm vector of permuted names
#' @param gs gene set
#' @param fract_min only cluster of size less or equal to this fraction of cell with not null feature values will be analysed
#' @description Calculate permutations

calc_gs_perm <- function(rll=NULL, perm=NULL, gs=NULL, fract_min=0.2){

  out <- setNames(numeric(length(rll)), names(rll))
  
  for(i in 1:length(rll)){
    
    #only if the current cluster_size is less or equal to fract_min*feature_size
    if(length(gs)*fract_min <= length(rll[[i]])){
      out[i] <- es(idx=which(perm[[i]] %in% gs), x=array(rll[[i]], dimnames=list(perm[[i]])))[, 1]
    }else{
      out[i] <- 0
    }
    
  }
  
  return(out)

}
