#' CNV Heatmap
#' @param scMuffinList scMuffinList object
#' @param file File name to save the figure as png file
#' @param width image width
#' @param height image height
#' @param units image units
#' @param res image resolution
#' @param genes gene labels to show as row annotation
#' @param genes.labels whether to show gene labels or just mark their position
#' @param mark.detected.cnv whether to show detected CNV regions
#' @param cluster_fontsize cluster fontsize
#' @param chrom_fontsize chromosome fontsize
#' @param legend_fontsize legend fontsize
#' @param genes.labels.fontsize gene labels fontsize
#' @param ... arguments passed to ComplexHeatmap::Heatmap
#' @importFrom ComplexHeatmap rowAnnotation anno_mark Heatmap draw
#' @importFrom grid gpar unit
#' @importFrom grDevices png
#' @export

heatmap_CNV <- function(scMuffinList = NULL, genes=NULL, genes.labels=FALSE, mark.detected.cnv=FALSE, file=NULL, width=180, height=180, units="mm", res=300, cluster_fontsize=8, chrom_fontsize=8, legend_fontsize=8, genes.labels.fontsize=8, ...) {
  
  #column_title_fontsize=8
  
  cnv_clustering <- scMuffinList$partitions[, "CNV"]
  names(cnv_clustering) <- rownames(scMuffinList$partitions)
  cnv <- scMuffinList$CNV$full$CNV
  
  #same cell ordering
  cnv <- cnv[, match(names(cnv_clustering), colnames(cnv))]
  
  ### add requested genes
  ha<-NULL
  if(!is.null(genes)){
    idx_genes <- setNames(vector("list", length(genes)), genes)
    for(i in 1:length(genes)){
      idx_genes[[i]] <- which(unlist(lapply(scMuffinList$CNV$full$regions2genes, function(x) any(x$symbol == genes[i]))))
      idx_genes[[i]] <- median(which(unlist(lapply(scMuffinList$CNV$full$regions2genes, function(x) any(x$symbol == genes[i])))))
      
    }
    idx_genes <- idx_genes[lengths(idx_genes)>0]
    if(length(idx_genes)>0){
      idx_cnv <- rep("", nrow(cnv))
      for(i in 1:length(idx_genes)){
        idx_cnv[idx_genes[[i]]] <- paste(idx_cnv[idx_genes[[i]]], names(idx_genes)[i], sep = ";")
      }
      idx_cnv <- gsub(";$", "", idx_cnv)
      idx_cnv <- gsub("^;", "", idx_cnv)
      idx_cnv <- factor(idx_cnv)
      
      if(genes.labels){
        ha <- rowAnnotation(link = anno_mark(at = which(idx_cnv != ""), labels = idx_cnv[idx_cnv != ""], labels_gp = gpar(fontsize = genes.labels.fontsize), padding = unit(1, "mm")))
      }else{
        #ha <- rowAnnotation(Gene = idx_cnv, annotation_legend_param=list(title_gp = gpar(fontsize = legend_fontsize), labels_gp = grid::gpar(fontsize = legend_fontsize, fontface = "italic")), show_annotation_name=FALSE)
        ha <- rowAnnotation(Gene = ifelse(idx_cnv=="", 0, 1), col=list(Gene=setNames(c("white", "black"), c(0, 1))))
      }
    }
  }
  
  if(mark.detected.cnv){
    idx_cnv <- rep(0, nrow(cnv))
    for(i in 1:length(scMuffinList$CNV$full$detected_cnv_regions)){
      for(j in 1:nrow(scMuffinList$CNV$full$detected_cnv_regions[[i]])){
        idx_start <- which(rownames(cnv) == scMuffinList$CNV$full$detected_cnv_regions[[i]]$start[j])
        idx_stop <- which(rownames(cnv) == scMuffinList$CNV$full$detected_cnv_regions[[i]]$stop[j])
        idx_cnv[idx_start:idx_stop] <- 1
      }
    }
    ha <- rowAnnotation(CNV = idx_cnv, col=list(CNV=setNames(c("white", "black"), c(0, 1))))
  }
  
  
  row_chr <- gsub("(chr[^_]+)_.+", "\\1", rownames(cnv))
  row_chr <- factor(row_chr, levels = unique(row_chr))
  row_chr <- split(row_chr, row_chr)
  ngenes_chr <- unlist(lapply(row_chr, length))
  row_splits <- factor(rep(names(ngenes_chr), ngenes_chr), levels = unique(names(ngenes_chr)))
  col_splits <- cnv_clustering
  clust_color <- rep("black", length(levels(cnv_clustering)))
  
  #ref_cluster <- cnv_clustering[names(cnv_clustering) == "reference"]
  ref_cluster <- scMuffinList$CNV$full$ref_cluster
  
  if(length(ref_cluster)>0){
    clust_color[levels(cnv_clustering)==ref_cluster] <- "red"
  }
  
  column_title_gp <- gpar(col = clust_color, fontsize = cluster_fontsize)
  row_title_gp <- gpar(fontsize = chrom_fontsize)
  heatmap_legend_param <- list(title_gp = gpar(fontsize = legend_fontsize), labels_gp = gpar(fontsize = legend_fontsize))
  
  max_abs <- max(abs(cnv), na.rm = T)
  
  if(!is.null(file)){
    png(file, width=width, height=height, units=units, res=res)
  }
  
  
  temp <- Heatmap(cnv, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, row_split = row_splits, row_title_rot=0, row_gap = unit(0, "mm"), column_split = col_splits, column_gap = unit(0, "mm"), border=TRUE, column_title_gp=column_title_gp, row_title_gp=row_title_gp, name="C", heatmap_legend_param=heatmap_legend_param, right_annotation = ha, ...)
  draw(temp)
  
  if(!is.null(file)){
    dev.off()
  }
  
  
}
