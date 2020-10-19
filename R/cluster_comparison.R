#' cluster_comparison
#' 
#' 
#' @export

cluster_comparison <- function(cl1, cl2, ctree1, ctree2){
	
	ctree1 <- as.hclust(ctree1)
	ctree2 <- as.hclust(ctree2)
	
	ans <- table(cl1, cl2)
	ans <- t(apply(ans, 1, function(x) x / colSums(ans)))
	ans <- ans[match(ctree1$labels, rownames(ans)), match(ctree2$labels, colnames(ans))]
	
	ans <- as.data.frame.ts(ans)
	
	return(ans)
}