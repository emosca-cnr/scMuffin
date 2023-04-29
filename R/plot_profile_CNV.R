#' CNV Profile
#' @param scMuffinList scMuffinList object
#' @param z.score whether to plot clustere median z-scores of CNV signal (TRUE) or cluster median CNV signal (FALSE).
#' @param cluster cluster to be plotted
#' @param file output file
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param image_format png or jpeg
#' @param cex.points cex of points
#' @param cex.lab cex of labs
#' @param cex.axis cex of axis
#' @param cex.main cex of main
#' @description Plot an heatmap of the CNV.
#' @details CNV Profile of every cluster
#' @import ComplexHeatmap grDevices grid
#' @export

plot_profile_CNV <- function(scMuffinList = NULL, cluster=0, z.score=TRUE, file=NULL, width=300, height=90, units="mm", res=300, image_format="png", cex.points=0.7, cex.lab=0.7, cex.axis=0.7, cex.main=0.7){
  
  cnv_clustering <- scMuffinList$partitions[, "CNV"]
  names(cnv_clustering) <- rownames(scMuffinList$partitions)
  cnv <- scMuffinList$CNV$full$CNV
  
  cnv <- cnv[, match(names(cnv_clustering), colnames(cnv))]
  
  row_chr <- gsub("(chr[^_]+)_.+", "\\1", rownames(cnv))
  row_chr <- factor(row_chr, levels = unique(row_chr))
  row_chr <- split(row_chr, row_chr)
  ngenes_chr <- unlist(lapply(row_chr, length))
  
  if(z.score){
    avg_cl_cnv <- t(apply(t(scale(t(cnv))), 1, function(i_row) tapply(i_row, cnv_clustering, median)))
  }else{
    avg_cl_cnv <- t(apply(cnv, 1, function(i_row) tapply(i_row, cnv_clustering, median)))
    avg_cnv <- apply(cnv, 1, median)
  }
  #sd_cl <- apply(avg_cl_cnv, 2, sd)
  sd_all <- sd(as.numeric(avg_cl_cnv))
  
  cl_lev <- levels(cnv_clustering)
  
  j <- which(cl_lev == cluster)
  #for(j in 1:length(cl_lev)){
    
    if(!is.null(file)){
      if(image_format == "jpeg"){
        file <- paste0(file, "_", cl_lev[j], ".jpg")
        jpeg(file, width=width, height=height, units=units, res=res)
      }
      if(image_format == "png"){
        file <- paste0(file, "_", cl_lev[j], ".png")
        png(file, width=width, height=height, units=units, res=res)
      }
    }
    
    #jpeg("prova.jpg", width=200, height=100, units="mm", res=300)
    
    par(mar=c(3, 3, 2,.1))
    par(mgp=c(2, 1, 0))
    
    ylim <- c(min(avg_cl_cnv[, j]), max(avg_cl_cnv[, j]))
    ylim <- c(min(avg_cl_cnv), max(avg_cl_cnv))
    plot(0, pch="", ylab="CNV signal", xlab="", xaxt="n", xlim=c(1, nrow(avg_cl_cnv)), ylim=ylim, xaxs="i", main=paste0("Cluster ", cl_lev[j]), cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main)
    xleft <- c(-10, cumsum(ngenes_chr)+0.5)
    xleft <- xleft[1:(length(xleft)-1)]
    xright <- cumsum(ngenes_chr)+0.5
    xright[length(xright)] <- xright[length(xright)]+10
    #rect(xleft, ylim[1]-1, xright, ylim[2]+2, col=adjustcolor(pals::alphabet(length(ngenes_chr)), 0.2))
    rect(xleft, ylim[1]-1, xright, ylim[2]+2, col="white")
    
    ### colors
    idx_neg <- avg_cl_cnv[, j] <= 0
    idx_pos <- avg_cl_cnv[, j] > 0
    md <- c(setNames(ggplot2::cut_interval(avg_cl_cnv[idx_neg, j], 5, dig.lab = 2), rownames(avg_cl_cnv)[idx_neg]), setNames(ggplot2::cut_interval(avg_cl_cnv[idx_pos, j], 5, dig.lab = 2), rownames(avg_cl_cnv)[idx_pos]))
    md <- md[match(rownames(avg_cl_cnv), names(md))]
    md <- factor(md, levels=levels(md))
    
    col <- rev(pals::brewer.rdylbu(10))[as.numeric(md)]
    
    chr_start <- c(1, cumsum(ngenes_chr)[-length(ngenes_chr)]+1)
    chr_end <- cumsum(ngenes_chr)
    for(k in 1:length(ngenes_chr)){
      #lines(chr_start[k]:chr_end[k], avg_cl_cnv[chr_start[k]:chr_end[k], j], col=col[(k %% 2) + 1])
      points(chr_start[k]:chr_end[k], avg_cl_cnv[chr_start[k]:chr_end[k], j], col=col[chr_start[k]:chr_end[k]], pch=16, cex=cex.points)
    }
    #abline(h=c(-sd_cl[j], sd_cl[j]), lty=2)
    abline(h=c(-sd_all, sd_all), lty=2)
    #lines(avg_cl_cnv[, j], col="red")
    #if(!z.score){
    #  lines(avg_cnv)
    #}
    
    abline(h=0, lty=2)
    axis(1, at = cumsum(ngenes_chr) - ngenes_chr/2, labels = names(ngenes_chr), cex.axis=cex.axis, las=2, cex.lab=cex.lab)
    
    if(!is.null(file)){
      dev.off()
    }
  #}

}
