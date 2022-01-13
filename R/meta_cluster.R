#' Returns meta-clusters obtained by using hierarchical clustering
#' @param overlap_mat Overlap matrix as obtained form overlap_matrix() function
#' @param n_step number used to divide hclust dendrogram height to obtain cut values to calculate silhouette
#' @param max_clust maximum number of clusters
#' @import cluster
#' @importFrom stats hclust as.dist cutree
#' @export

meta_cluster <- function(overlap_mat, n_step = 11, max_clust=10) {
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
  plot(sil$step, sil$silhouette, type="b", main="Silhouette",
       xlab="steps", ylab="Average Silhouette")
  abline(v=sil$step[which(sil$silhouette == max(sil$silhouette, 
                                                      na.rm = T))[1]], 
         col="red",
         lwd=3, lty=2)
  
  clusters_hc <- stats::cutree(hc_cl, h = sil$step[which(sil$silhouette == max(sil$silhouette, na.rm = T))[1]])
  n_clust <- max(clusters_hc)
  if(n_clust>max_clust){
  	clusters_hc <- cutree(hc_cl, k = max_clust)
  }
  clusters <- data.frame(cluster= names(clusters_hc), meta_cl = clusters_hc, stringsAsFactors = F)
  ans <- list(dissimilarity = diss, hclust_out = hc_cl, clusters = clusters)
  
  return(ans)
}




