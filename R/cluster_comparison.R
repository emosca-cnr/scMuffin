#' cluster_comparison
#' @param cl1 x
#' @param cl2 x
#' @param ctree1 x
#' @param ctree2 x
#' @importFrom stats as.hclust
#' @description Cluster comparison
#' 
#' @export

cluster_comparison <- function(cl1, cl2, ctree1, ctree2){
	
	ctree1 <- stats::as.hclust(ctree1)
	ctree2 <- stats::as.hclust(ctree2)
	
	ans <- table(cl1, cl2)
	ans_freq <- t(apply(ans, 1, function(x) x / colSums(ans)))
	
	ans <- ans[match(ctree1$labels, rownames(ans)), match(ctree2$labels, colnames(ans))]
	ans_freq <- ans_freq[match(ctree1$labels, rownames(ans_freq)), match(ctree2$labels, colnames(ans_freq))]
	
	ans <- as.data.frame.ts(ans)
	ans_freq <- as.data.frame.ts(ans_freq)
	
	return(list(occ=ans, freq=ans_freq))
}