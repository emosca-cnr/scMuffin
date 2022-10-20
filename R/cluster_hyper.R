#' Calculate cluster enrichment by hypergeometric test
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param feature_values matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param fdr fdr threshold, default at 0.05
#' @param mc.cores number of cores
#' @importFrom parallel mclapply
#' @description Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @export
#' @author Ettore Mosca

cluster_hyper <- function(feature_values, cell_clusters, fdr=0.05, top=2, mc.cores=1){
  
  
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
      ora_res[[i]][[j]]  <- parallel::mclapply(gsl, function(x) ora(feature_val_i[[j]], universe[!universe %in% feature_val_i[[j]]], gsl = x, p_adj_method = "fdr"), mc.cores = mc.cores)
      
    }
    
  }
  
  return(ora_res)
}
  