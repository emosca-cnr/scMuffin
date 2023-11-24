#' Enrichment Score
#'
#' @param idx vector of indices a subset of elements of x
#' @param x named vector, ranked list
#' @return data.frame with as, enrichmet score; tags, leading edge size; tags_perc, leading edge size percent over gene set; list_top, rank of the ES; list_top_perc, rank of the ES percent over full ranked list; lead_edge, gene names of the leading edge.
#' @export
#'
es <- function(idx=NULL, x=NULL){
  
  #idx: indexes of a subset of elements of x
  #x: array of elements
  ## ES score
  #Phit(S,i) <- SUM_{i in idx, j=1..i}( r_j / Nr)
  #r: score
  #Nr: total score in x[idx]
  
  N <- length(x)
  Nh <- length(idx)
  hits <-  rep(0, length(x))
  hits[idx] <- x[idx]
  misses <- rep(1, length(x))
  misses[idx] <- 0
  hits.cumsum <- cumsum(abs(hits))
  Nr <- sum(abs(hits)) #equal to sum(abs(x[idx]))
  
  if(Nr==0){   #nothing to do
    
    es <- 0
    
    deviation <- 0
    tags <- 0
    tags_perc <- 0
    list_top <- 0
    list_top_perc <- 0
    lead_edge <- 0
    lead_edge_subset <- 0
    
  }else{
    
    misses.cumsum <- cumsum(misses)
    Phits <- hits.cumsum / Nr
    Pmiss <- misses.cumsum / (N - Nh)
    deviation <- Phits - Pmiss
    wm <- which.max(abs(deviation))
    
    #es <- deviation[which.max(abs(deviation))]
    es <- deviation[wm]
    
    
    if(es >=0){
      
      tags <- sum(idx <= wm)
      list_top <- wm
      lead_edge_subset <- paste0(names(x)[intersect(1:wm, idx)], collapse = ";")
      
    }else{
      
      tags <- sum(idx >= wm)
      list_top <- N - (wm-1)
      lead_edge_subset <- paste0(names(x)[intersect(wm:N, idx)], collapse = ";")
      
    }
    
    tags_perc <- tags / length(idx)
    list_top_perc <- list_top / length(x)
    
  }
  
  lead_edge <- tags_perc * (1-list_top_perc) * N / (N - Nh)
 
  return(data.frame(es=es, tags=tags, tags_perc=tags_perc, list_top=list_top, list_top_perc=list_top_perc, lead_edge=lead_edge, lead_edge_subset=lead_edge_subset, stringsAsFactors = F))
  
}
