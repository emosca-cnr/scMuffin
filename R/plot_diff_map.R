#' Plot diffusion map colored by other features
#' @param scMuffinList scMuffinList object
#' @param columns columns of scMuffinList$diffusion_map_pseudo_t$summary to be used for x and y
#' @param col_data named vector with values that will be translated into colors. Names are cell names.
#' @param file File name to save the figure as png file
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param min_cells minimum number of cells in which the feature must have a non-zero value
#' @param ... further arguments to plot()
#' @param scale_feature logical, whether to scale col_data
#' @param adj_outliers logical, whether to adjust the group.by scores, removing outliers
#' @description Produce a scatter plot where cells are placed according to the results of diffusion map analysis and colored by values given in col_data.
#' @import grDevices 
#' @export

plot_diff_map <- function(scMuffinList = NULL, columns=c(1:2), col_data=NULL, file=NULL, width=200, height=200, units="mm", res=300, adj_outliers=TRUE, scale_feature=FALSE, min_cells=50, ...){
  
  if(length(scMuffinList$diffusion_map_pseudo_t)==0){
    stop("scMuffinList$diffusion_map_pseudo_t does not exist\n")
  }
  
  if(!is.null(col_data)){
    if(!all(rownames(scMuffinList$diffusion_map_pseudo_t$summary) %in% names(col_data))){
      stop("No all rownames(scMuffinList$diffusion_map_pseudo_t$siummary) were found in names(col_data).\n")
    }else{
      col_data <- col_data[match(rownames(scMuffinList$diffusion_map_pseudo_t$summary), names(col_data))]
    }
  }else{
    stop("col_data is required.\n")
  }
  
  
  if(!is.null(file)){
    png(file, width=width, height=height, units=units, res=res)
  }
  
  #jpeg("prova.jpg", width=200, height=100, units="mm", res=300)
  
  
  if(is.numeric(col_data)){
    
    if(adj_outliers){
      col_data <- adj_outliers_col(col_data)
      cat("Adjusting data for outliers\n")
    }
    
    if(scale_feature){
      col_data <- setNames(as.numeric(scale(col_data)), names(col_data))
      cat("Scaling data\n")
      
    }
    
    col_data[is.na(col_data)] <- 0
    
    if(sum(abs(col_data) > 0) >= min_cells){
      
      if(any(col_data < 0)){
        
        #cut the feature id separately for <=0 and >0 in a total of 10 intervals
        #reorder cells according to feature data
        idx_neg <- col_data <= 0
        idx_pos <- col_data > 0
        md <- c(setNames(ggplot2::cut_interval(col_data[idx_neg], 5, dig.lab = 2), names(col_data)[idx_neg]), setNames(ggplot2::cut_interval(col_data[col_data>0], 5, dig.lab = 2), names(col_data)[idx_pos]))
        md <- md[match(names(col_data), names(md))]
        md <- factor(md, levels=rev(levels(md)))
        pal <- pals::brewer.rdylbu(10)
        cols <- pal[as.numeric(md)]
        
      }else{
        
        #if only positive values
        md <- ggplot2::cut_interval(col_data, 5, dig.lab = 2)
        md <- factor(md, levels=rev(levels(md)))
        pal <- rev(pals::brewer.ylorrd(5))
        cols <- pal[as.numeric(md)]
        
      }
    }
    
  }else{        #if not numeric...
    
    md <- as.factor(col_data)
    n_colors <- length(levels(md))
    pal <- pals::alphabet(n_colors)
    cols <- pal[as.numeric(as.factor(col_data))]
    
  }
  
  layout(matrix(1:2, nrow = 1), widths = c(0.8, 0.2))
  par(mar=c(3, 3, 2, .1))
  par(mgp=c(2, 1, 0))
  
  xx <- scMuffinList$diffusion_map_pseudo_t$summary[, 1]
  yy <- scMuffinList$diffusion_map_pseudo_t$summary[, 2]
  
  plot(xx, yy, col=cols, pch=16, xlab=columns[1], ylab=columns[2], ...)
  
  par(mar=c(3, .1, 2, .1))
  plot.new()
  legend("center", legend = levels(md), pch=16, col=pal, cex = 0.5)
  
  
  if(!is.null(file)){
    dev.off()
  }
  
}
