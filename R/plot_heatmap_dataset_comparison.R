#' plot_heatmap_dataset_comparison
#' @param dataset_cmp_list list resulting from inter_dataset_comparison()
#' @param type "score" or "significance" for, respectively, cluster median marker set score or -log10(FDRq) of CSEA.
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image res
#' @param outfile File name to save the figure as png file
#' @param show.gs.source whether to show or not gene set source
#' @param ... arguments passed to ComplexHeatmap::Heatmap
#' @export
#' @import ComplexHeatmap pals
#' @importFrom circlize colorRamp2

plot_heatmap_dataset_comparison <- function(dataset_cmp_list=NULL, type="score", show.gs.source=FALSE, outfile=NULL, width=200, height=200, units="mm", res=300, ...){
  
  type <- match.arg(type, c("score", "significance"))
  if(!is.null(outfile)){
    png(outfile, width = width, height = height, res=res, units = units)
  }
  
  if(type=="score"){
    M <- dataset_cmp_list$score_matrix
    pal <- rev(pals::brewer.rdylbu(11))
    name <- "score"
  }else{
    M <- -log10(dataset_cmp_list$significance_matrix)
    pal <- pals::brewer.purples(11)
    name <- "sig"
  }
  M <- as.matrix(M)
  X <- adj_outliers_col(as.numeric(M))
  
  if(any(X<0)){
    X_extreme <- max(abs(min(X)), max(X))
    X_extreme <- c(-X_extreme, X_extreme)
  }else{
    X_extreme <- max(X)
    X_extreme <- c(0, X_extreme)
  }
  
  col_fun <- circlize::colorRamp2(seq(X_extreme[1], X_extreme[2], length.out = 11), pal)
  
  ra <- factor(gsub("^([^_]+)_.+", "\\1", rownames(M)))
  ra <- rowAnnotation(clusters=ra, col=list(clusters=setNames(pals::brewer.accent(length(levels(ra)))[as.numeric(ra)], ra)))
  ca <- NULL
  if(show.gs.source){
    ca <- factor(gsub("^([^_]+)_.+", "\\1", colnames(M)))
    ca <- columnAnnotation(markers=ca, col=list(markers=setNames(pals::brewer.pastel1(length(levels(ca)))[as.numeric(ca)], ca)))
  }
  
  hm <- Heatmap(as.matrix(M), col=col_fun, right_annotation = ra, top_annotation = ca, name = name, ...)
  draw(hm)
  
  if(!is.null(outfile)){
    dev.off()
  }
  #return(hm)
  
}