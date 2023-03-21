#' Cell cluster statistics
#' @param scMuffinList scMuffinList object
#' @param partition_id partition id
#' @param feature_id feature id
#' @param mean_f function to use for average
#' @param var_f function to use for variability
#' @param na.rm whether to remove na or not
#' @importFrom stats sd setNames
#' @description Calculate mean and variance for partition id and feature id
#' @return scMuffinList with data.frames "mean" and "var" added to "scMuffinList$cluster_data[[partition_id]]"
#' @export

cluster_stats <- function(scMuffinList=NULL, partition_id=NULL, feature_id=NULL, mean_f=mean, var_f=sd, na.rm=TRUE){
  
  X <- as.matrix(scMuffinList[[feature_id]]$summary)
  X_clusters <- setNames(scMuffinList$partitions[, partition_id], rownames(scMuffinList$partitions))
  
  X_clusters <- X_clusters[match(rownames(X), names(X_clusters))]
  
  if(is.numeric(X)){
    X_mean <- apply(X, 2, function(i_col) tapply(i_col, X_clusters, FUN = mean_f, na.rm=na.rm))
    X_var <- apply(X, 2, function(i_col) tapply(i_col, X_clusters, FUN = var_f, na.rm=na.rm))
  }else{
    stop("scMuffinList[[feature_id]]$summary is not numeric\n")
  }
  
  if(length(scMuffinList$cluster_data[[partition_id]]$mean)==0){
    
    scMuffinList$cluster_data[[partition_id]]$mean <- X_mean
    scMuffinList$cluster_data[[partition_id]]$var <- X_var
    
  }else{
    
    idx_shared <- intersect(colnames(scMuffinList$cluster_data[[partition_id]]$mean), colnames(X))
    
    if(length(idx_shared)>0){
      scMuffinList$cluster_data[[partition_id]]$mean[, idx_shared] <- NULL
      scMuffinList$cluster_data[[partition_id]]$var[, idx_shared] <- NULL
    }
    
    scMuffinList$cluster_data[[partition_id]]$mean <- merge(scMuffinList$cluster_data[[partition_id]]$mean, X_mean, by=0)
    rownames(scMuffinList$cluster_data[[partition_id]]$mean) <- scMuffinList$cluster_data[[partition_id]]$mean$Row.names
    scMuffinList$cluster_data[[partition_id]]$mean$Row.names <- NULL
    
    scMuffinList$cluster_data[[partition_id]]$var <- merge(scMuffinList$cluster_data[[partition_id]]$var, X_var, by=0)
    rownames(scMuffinList$cluster_data[[partition_id]]$var) <- scMuffinList$cluster_data[[partition_id]]$var$Row.names
    scMuffinList$cluster_data[[partition_id]]$var$Row.names <- NULL
    
    
  }
  
  return(scMuffinList)
  
}