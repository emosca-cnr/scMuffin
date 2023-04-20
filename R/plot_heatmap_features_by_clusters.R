#' Plot an heatmap of feature values in cell clusters
#' @param scMuffinList scMuffinList object
#' @param significance_matrix optional significance matrix (clusters-by-features) of the same size of the data specified by means of feature_source
#' @param feature_source It can be a "mean", "gss" or a numeric matrix (clusters-by-features). If "mean", the data.frame with average feature values among clusters will be used (default); if "gene_set_scoring", the average gene set values among clusters will be used.
#' @param sig_threshold significance threshold
#' @param file output file
#' @param remove_null_features whether to remove null features
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param partition_id identifier of the partition to be considered
#' @param image_format png or jpeg
#' @param scale whether to scale the features
#' @param pal color palette. Default to rev(pals::brewer.rdylbu(10)) (negative values) or pals::brewer.ylorrd(5)) (positive values)
#' @param ... further arguments to ComplexHeatmap::Heatmap
#' @export
#' @import ComplexHeatmap grDevices
#' @importFrom utils write.table
#' @importFrom pals brewer.rdylbu brewer.ylorrd
#' @importFrom circlize colorRamp2

plot_heatmap_features_by_clusters <- function(scMuffinList=NULL, feature_source=NULL, partition_id=NULL, significance_matrix=NULL, sig_threshold=0.05, file=NULL, remove_null_features=FALSE, width=180, height=180, units="mm", res=300, scale=TRUE, pal=NULL, ...){
  
  
  if(!is.null(significance_matrix)){
    cell_fun_asterisk <- function(j, i, x, y, w, h, fill) {
      if(sig_mat[i, j] < sig_threshold) {
        grid.text("*", x, y)
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

  if(remove_null_features){
    X[is.na(X)] <- 0
    X <- X[rowSums(abs(X))>0, ]
    print(dim(X))
  }
  
  if(is.null(pal)){
    X_abs_max <- max(abs(c(min(X, na.rm = TRUE), max(X, na.rm = TRUE))))
    if(any(X<0)){
      pal <- rev(pals::brewer.rdylbu(5))
      pal <- colorRamp2(c(-X_abs_max, 0, X_abs_max), c(pal[1], pal[3], pal[5]))
    }else{
      pal <- pals::brewer.ylorrd(5)
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
    }
    
    h_tot_go <- ComplexHeatmap::Heatmap(X, show_row_names = T, cell_fun = cell_fun_asterisk, col=pal, ...)
    draw(h_tot_go, heatmap_legend_side = "left")
    
    if(!is.null(file)){
      dev.off()
    }
    

  }else{
    message("Are all rows empty?\n")
  }
  
}


