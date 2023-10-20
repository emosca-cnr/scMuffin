#' Boxplot with points
#' @param x values
#' @param f factor to split the values
#' @param col color palette. Equal in length to levels(f)
#' @param amount amount of jitter
#' @param adj.col color adjustment
#' @param pch point pch
#' @param cex cex
#' @param ylim ylim
#' @param file output file
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param image_format png or jpeg
#' @param ... further argument to boxplot
#' @importFrom pals alphabet
#' @importFrom grDevices png jpeg
#' @export

boxplot_points <- function(x=NULL, f=NULL, col=NULL, amount=0.2, adj.col=1, pch=16, cex=0.6, ylim=NULL, file=NULL, width=180, height=180, units="mm", res=300, image_format="png", ...){
  
  f <- as.factor(f)
  f_lev <- levels(f)
  
  if(is.null(col)){
    col <- alphabet(length(f_lev))
  }
  
  if(is.null(ylim)){
    ylim=c(min(x), max(x))
  }
  
  
  if(!is.null(file)){
    if(image_format == "jpeg"){
      jpeg(file, width=width, height=height, units=units, res=res)
    }
    if(image_format == "png"){
      png(file, width=width, height=height, units=units, res=res)
    }
  }
  
  
  boxplot(x ~ f, outline=F, col=NA, border=col, ylim=ylim, ...)
  
  for(i in 1:length(f_lev)){
    xx <- jitter(rep(i, sum(f==f_lev[i])), amount=amount)
    yy <- x[f==f_lev[i]]
    points(xx, yy, col=adjustcolor(col[i], adj.col), pch=pch, cex=cex)
  }
  
  if(!is.null(file)){
    dev.off()
  }
}