#' Plot an heatmap of feature values in cell clusters
#' @param scMuffinList scMuffinList object
#' @param significance_matrix optional significance matrix (clusters-by-features) of the same size of the data specified by means of feature_source
#' @param feature_source It can be a "mean", "gss" or a numeric matrix (clusters-by-features). If "mean", the data.frame with average feature values among clusters will be used (default); if "gene_set_scoring", the average gene set values among clusters will be used.
#' @param sig_threshold significance threshold
#' @param file File name to save the figure as png file
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param partition_id identifier of the partition to be considered
#' @param scale whether to scale the features
#' @param pal color palette. Default to rev(pals::brewer.rdylbu(10)) (negative values) or pals::brewer.ylorrd(5)) (positive values)
#' @param na_col color for NA values
#' @param X_abs_max maximum absolute value permitted, useful to avoid the effect of outliers over colors
#' @param ... further arguments to ComplexHeatmap::Heatmap
#' @export
#' @import ComplexHeatmap grDevices pals
#' @importFrom circlize colorRamp2

plot_heatmap_features_by_clusters <- function(scMuffinList=NULL, feature_source=NULL, partition_id=NULL, significance_matrix=NULL, sig_threshold=0.05, file=NULL, width=180, height=180, units="mm", res=300, scale=FALSE, pal=NULL, na_col="black", X_abs_max=NULL, ...){
  
  
  if(!is.null(significance_matrix)){
    cell_fun_asterisk <- function(j, i, x, y, w, h, fill) {
      if(!is.na(sig_mat[i, j])){
        if(sig_mat[i, j] < sig_threshold) {
          grid.text("*", x, y)
        }
      }
    }
  }else{
    cell_fun_asterisk <- NULL
  }
  
  
  if(is.matrix(feature_source)){
    X <- t(feature_source)
  }else{
    
    if(length(scMuffinList$cluster_data)==0){
      stop("If feature_source is not a matrix, then scMuffinList$cluster_data must be defined.\n")
    }
    
    partition_id <- match.arg(partition_id, names(scMuffinList$cluster_data))
    
    if(is.null(feature_source)){
      X <-  t(scMuffinList$cluster_data[[partition_id]]$mean)
    }else{
      if(feature_source == "mean"){
        X <-  t(scMuffinList$cluster_data[[partition_id]]$mean)
      }
      if(feature_source == "gss"){
        X <-  t(as.matrix(scMuffinList$cluster_data[[partition_id]]$gene_set_scoring$summary))
      }
    }
  }
  
  #Detect and remove null features
  null_rows <- which(rowSums(abs(X), na.rm = T)==0)
  if(length(null_rows) > 0){
    cat("Detected features with all null values: these will be removed...\n")
    if(length(null_rows) == nrow(X)){
      stop("All features are null.\n")
    }
    X <- X[-null_rows, , drop=FALSE]
  }
  
  # if(remove_null_features){
  #   X[is.na(X)] <- 0
  #   X <- X[rowSums(abs(X))>0, ]
  #   print(dim(X))
  # }
  
  if(is.null(pal)){
    if(is.null(X_abs_max)){
      X_abs_max <- max(abs(c(min(X, na.rm = TRUE), max(X, na.rm = TRUE))))
    }else{
      X[X > X_abs_max] <- X_abs_max
      X[X < -X_abs_max] <- -X_abs_max
    }
    if(any(X<0)){
      pal <- rev(brewer.rdylbu(5))
      pal <- colorRamp2(c(-X_abs_max, 0, X_abs_max), c(pal[1], pal[3], pal[5]))
    }else{
      pal <- brewer.ylorrd(5)
      pal <- colorRamp2(c(0, X_abs_max), c(pal[1], pal[5]))
    }
  }
  
  
  if(nrow(X)>0){
    
    if(!is.null(significance_matrix)){
      sig_mat <- t(significance_matrix)
      sig_mat <- sig_mat[match(rownames(X), rownames(sig_mat)), match(colnames(X), colnames(sig_mat))]
      # if(onlyUp){
      #   sig_mat[!X>0] <- 1
      # }
    }
    
    if(!is.null(file)){
      png(file, width=width, height=height, units=units, res=res)
    }
    
    if(scale){
      X <- t(scale(t(X)))
      cat("Distribution of z-scores, set X_abs_max to control the colors of the maximum absolute values:\n")
      print(summary(as.numeric(unlist(X))))
    }
    
    h_tot_go <- Heatmap(X, show_row_names = T, cell_fun = cell_fun_asterisk, col=pal, na_col=na_col, ...)
    draw(h_tot_go, heatmap_legend_side = "left")
    
    if(!is.null(file)){
      dev.off()
    }
    
    
  }else{
    message("Are all rows empty?\n")
  }
  
}


