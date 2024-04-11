#' plot_heatmap_dataset_comparison
#' @param dataset_cmp_list list resulting from inter_dataset_comparison()
#' @param type "score" or "significance" for, respectively, cluster median marker set score or -log10(FDRq) of CSEA.
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image res
#' @param outfile File name to save the figure as png file
#' @param show.gs.source whether to show or not gene set source
#' @param na_col color for NA values
#' @param ... arguments passed to ComplexHeatmap::Heatmap
#' @importFrom ComplexHeatmap rowAnnotation columnAnnotation Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom pals brewer.pastel1 brewer.accent brewer.rdylbu brewer.purples
#' @export
#'  
plot_heatmap_dataset_comparison <- function(dataset_cmp_list=NULL, type="score", show.gs.source=FALSE, outfile=NULL, width=200, height=200, units="mm", res=300, na_col="black", ...){
  
  type <- match.arg(type, c("score", "significance"))
  if(!is.null(outfile)){
    png(outfile, width = width, height = height, res=res, units = units)
  }
  
  if(type=="score"){
    M <- dataset_cmp_list$score_matrix
    pal <- rev(brewer.rdylbu(11))
    name <- "score"
  }else{
    M <- -log10(dataset_cmp_list$significance_matrix)
    pal <- brewer.purples(11)
    name <- "sig"
  }
  M <- as.matrix(M)
  X <- adj_outliers_col(as.numeric(M))
  
  if(any(X<0, na.rm = T)){
    X_extreme <- max(abs(min(X, na.rm = T)), max(X, na.rm = T), na.rm = T)
    X_extreme <- c(-X_extreme, X_extreme)
  }else{
    X_extreme <- max(X, na.rm = T)
    X_extreme <- c(0, X_extreme)
  }
  
  col_fun <- colorRamp2(seq(X_extreme[1], X_extreme[2], length.out = 11), pal)
  
  ra <- factor(gsub("^([^_]+)_.+", "\\1", rownames(M)))
  ra <- rowAnnotation(clusters=ra, col=list(clusters=setNames(brewer.accent(length(levels(ra)))[as.numeric(ra)], ra)))
  ca <- NULL
  if(show.gs.source){
    ca <- factor(gsub("^([^_]+)_.+", "\\1", colnames(M)))
    ca <- columnAnnotation(markers=ca, col=list(markers=setNames(brewer.pastel1(length(levels(ca)))[as.numeric(ca)], ca)))
  }
  
  hm <- Heatmap(as.matrix(M), col=col_fun, right_annotation = ra, top_annotation = ca, name = name, na_col=na_col, ...)
  draw(hm)
  
  if(!is.null(outfile)){
    dev.off()
  }
  #return(hm)
  
}