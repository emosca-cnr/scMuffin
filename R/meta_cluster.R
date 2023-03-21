#' Analysis of the various partitions and definition of meta-clusters
#' @param scMuffinList scMuffinList object with at least two partitions
#' @param n_step number used to divide hclust dendrogram height to obtain cut values to calculate silhouette
#' @param max_clust maximum number of clusters
#' @param do_plot whether to do the silohouette plot or not
#' @import cluster
#' @importFrom stats hclust as.dist cutree
#' @description Calculate the overlap matrix between all-pars of clusters of existing partitions.
#' @return
#' scMuffinList with:
#' \itemize{
#'   \item{scMuffinList$cluster_comparison$overlap_matrix, the overlap matrix}
#'   \item{scMuffinList$cluster_comparison$meta_clusters, a list with information on meta-clusters:}
#'   \itemize{
#'     \item{dissimilarity, dissimilarity matrix;}
#'     \item{hclust_out, stats::hclust output;}
#'     \item{clusters, a data.frame for mapping  clusters and meta-clusters;}
#'     \item{cells, a data.frame for mapping cells and meta-clusters;}
#'     \item{occurrence, a list of meta-clusters, in which each element is a data.frame that show the occurrence of meta-cluster cells across clusters.}
#'     }
#' }
#' @export

meta_cluster <- function(scMuffinList=NULL, n_step = 11, max_clust=10, do_plot=FALSE) {
  

  scMuffinList$partitions$all <- as.factor(as.character(apply(scMuffinList$partitions[, colnames(scMuffinList$partitions) != "all"], 1, function(i_row) paste0(colnames(scMuffinList$partitions), i_row, collapse="__"))))
  

  scMuffinList <- overlap_matrix(scMuffinList = scMuffinList)
  overlap_mat <- scMuffinList$cluster_comparison$overlap_matrix
  
  diss <- max(overlap_mat) - overlap_mat
	hc_cl <- stats::hclust(stats::as.dist(diss))
	
	val <- seq(min(hc_cl[["height"]]), max(hc_cl[["height"]]), length.out =n_step)
	sil <- data.frame(step=val, silhouette= rep(0, length(val)), stringsAsFactors = F)
	for (i in 1:length(val)) {
		tmp <- cluster::silhouette(cutree(hc_cl, h=val[i]), diss)
		if(unique(as.vector(unlist(is.na(tmp))))) {
			sil[i,2] <- NA
		} else {
			sil[i,2] <- round(mean(tmp[,3]), digits = 2)
		}
	}
	
	if(do_plot){
		jpeg("silhouette.jpg", width = 90, height = 90, res=300, units="mm")
		plot(sil$step, sil$silhouette, type="b", main="Silhouette", xlab="steps", ylab="Average Silhouette")
		abline(v=sil$step[which(sil$silhouette == max(sil$silhouette, na.rm = T))[1]], col="red", lwd=3, lty=2)
		dev.off()
	}
	
	clusters_hc <- stats::cutree(hc_cl, h = sil$step[which(sil$silhouette == max(sil$silhouette, na.rm = T))[1]])
	n_clust <- max(clusters_hc)
	if(n_clust>max_clust){
		clusters_hc <- cutree(hc_cl, k = max_clust)
	}
	clusters <- data.frame(cluster= names(clusters_hc), meta_cl = clusters_hc, stringsAsFactors = F)
	
	
	scMuffinList$cluster_comparison$meta_clusters <- list(dissimilarity = diss, hclust_out = hc_cl, clusters = clusters)
	
	meta_clusters <- scMuffinList$cluster_comparison$meta_clusters
	cl_list <- scMuffinList$cluster_comparison$cluster_list
	
	meta_clusters$clusters$type <- gsub("_[^_]+$", "", meta_clusters$clusters$cluster)
	
	ans <- vector("list", nrow(meta_clusters$clusters))
	
	for(i in 1:nrow(meta_clusters$clusters)){
	  
	  idx_cl_list <- meta_clusters$clusters$type[i]
	  ans[[i]] <- data.frame(meta_cl=meta_clusters$clusters$meta_cl[i], type=meta_clusters$clusters$type[i], cluster_id=meta_clusters$clusters$cluster[i], cell_id=names(cl_list[[idx_cl_list]])[cl_list[[idx_cl_list]] == meta_clusters$clusters$cluster[i]], stringsAsFactors = F)
	  
	}
	
	ans <- Reduce(rbind, ans)
	scMuffinList$cluster_comparison$meta_clusters$cells <- unique(ans[, c("meta_cl", "cell_id")])
	
	temp <- split(ans, ans$meta_cl)
	temp <- lapply(temp, function(x) table(x$cell_id, x$cluster_id))
	scMuffinList$cluster_comparison$meta_clusters$occurrence <- temp

	return(scMuffinList)
	
}




