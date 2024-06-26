#' Calculate cluster enrichment by hypergeometric test
#' 
#' @param feature_values matrix, cells-by-features
#' @param cell_clusters array with cell values named by their cluster ID
#' @importFrom parallel mclapply
#' @description Calculate cluster enrichment in each of the categorical values of a feature, by means of hypergeometric test
#' @return A list of data.frames, one for each categorical value. See [ora()]
#' @export
#' @importFrom stats setNames

cluster_hyper <- function(feature_values=NULL, cell_clusters=NULL){
  
  
  #gene sets are clusters
  gsl <- list(clusters=lapply(split(cell_clusters, cell_clusters), function(x) names(x)))
  
  universe <- rownames(feature_values)
  
  ora_res <- vector("list", length = ncol(feature_values))
  names(ora_res) <- colnames(feature_values)
  for(i in 1:ncol(feature_values)){
    
    cat("processing", colnames(feature_values)[i], "...\n")
    feature_val_i <- lapply(split(setNames(feature_values[, i], rownames(feature_values)), feature_values[, i]), function(x) names(x))
    ora_res[[i]] <- vector("list", length(feature_val_i))
    names(ora_res[[i]]) <- names(feature_val_i)
    
    for(j in 1:length(feature_val_i)){
      ora_res[[i]][[j]]  <- ora(feature_val_i[[j]], universe[!universe %in% feature_val_i[[j]]], gsl = gsl$clusters, p_adj_method = "fdr")
    }
    
  }
  
  return(ora_res)
}
