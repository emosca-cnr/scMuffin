#' Co-occurrence of cells between two partitions
#' @description Calculate, as Jaccard index, the co-occurrence of cells in all the pairs of elements composed by an element from p1 and an element from p2.
#' @param p1 a partition among scMuffinList$partitions, e.g. scML_demo$partitions$global_exp
#' @param p2 a partition among scMuffinList$partitions, e.g. scML_demo$partitions$CNV
#' @return matrix of Jaccard indeces between all levels of p1 and all levels of p2
#' @export

cel_coocc_partitions <- function(p1=NULL, p2=NULL){
	
	stopifnot(length(p1) == length(p2))
	
	p1_by_p2 <- table(p1, p2)
	p1_size <- matrix(rep(rowSums(p1_by_p2), ncol(p1_by_p2)), ncol = ncol(p1_by_p2))
	p2_size <- matrix(rep(colSums(p1_by_p2), nrow(p1_by_p2)), nrow = nrow(p1_by_p2), byrow = T)
	p1_by_p2 <- p1_by_p2 / (p1_size + p2_size - p1_by_p2)
	p1_by_p2[is.nan(p1_by_p2)] <- 0
	
	return(p1_by_p2)
}
