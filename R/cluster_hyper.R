#' Calculate cluster enrichment by hypergeometric test
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features_by_cells matrix, features by cells matrix
#' @param cell_clusters array with cell values named by their cluster ID
#' @param n_top numeric, number of features to be shown with a difference color and representing the most significative features according to t-test
#' @param fdr fdr threshold, default at 0.05
#' @param dir_out string, output directory
#' @importFrom parallel mclapply
#' @export
#' @author Ettore Mosca

cluster_chisq <- function(feature_values, cell_clusters, fdr=0.05, top=2){
  
  
  #gene sets are clusters
  gsl <- list(clusters=lapply(split(cell_clusters, cell_clusters), function(x) names(x)))
  
  universe <- rownames(feature_values)
  
  #deg list are feature values
  for(i in 1:ncol(feature_values)){
    
    feature_val_i <- lapply(split(setNames(feature_values[, i], rownames(feature_values)), feature_values[, i]), function(x) names(x))
    ora_res <- vector("list", length(feature_val_i))
    names(ora_res) <- names(feature_val_i)
    for(j in 1:length(feature_val_i)){
      
      ora_res[[j]]  <- parallel::mclapply(gsl, function(x) ora(feature_val_i[[j]], universe[!universe %in% feature_val_i[[j]]], gsl = x, p_adj_method = "fdr"), mc.cores = mc.cores)
      
    }
    
  }
}
  