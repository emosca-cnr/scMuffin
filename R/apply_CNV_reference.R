#' apply_CNV_reference
#' @param cnv the output of calculate_CNV
#' @param cnv_clustering the output of cluster_by_features(..., cnv=TRUE, ...)
#' @param reference name of the reference column
#' @param eps threshold for the min_max method
#' @param method mean: subtract the average profile of the reference cluster to every cell; min_max see Tirosh et al.
#' @export

apply_CNV_reference <- function(cnv=NULL, cnv_clustering=NULL, reference="reference", eps=0.2, method=c("mean", "min_max")){
	
	method <- method[1]
	#cnv_clustering$clusters <- cnv_clustering$clusters[order(cnv_clustering$clusters)]
	#cnv <- cnv[, match(names(cnv_clustering$clusters), colnames(cnv))]
	
	#find the cluster where the reference appears
	if(is.null(reference) | !any(colnames(cnv) %in% reference)){
		stop("can't find the refence\n")
	}
	
	ref_cluster <- NULL
	if(!is.null(cnv_clustering)){
		
		ref_cluster <- cnv_clustering$clusters[names(cnv_clustering$clusters) == reference] #cluster in which the reference occurs
		idx_ref_cells <- which(colnames(cnv) %in% names(cnv_clustering$clusters)[cnv_clustering$clusters==ref_cluster])
		
		cat("Reference cluster:", as.character(ref_cluster), "\n")
		
	}
	
	if(length(reference) > 1){
		
		idx_ref_cells <- which(colnames(cnv) %in% reference)
		
		cat("Reference cells:", length(idx_ref_cells), "\n")
		
	}
	
	
	cat("Subtracting reference cells average from CNV profiles...\n")
	
	#update the CNV Matrix, subtracting the average of the reference cluster from CNV profiles
	if(method=="mean"){
		
		ref_cluster_avg <- rowMeans(cnv[, idx_ref_cells])
		cnv <- apply(cnv, 2, function(x) x-ref_cluster_avg)
		
	}
	if(method=="min_max"){
		
		#ref_cluster_min <- rowMin(cnv[, idx_ref_cells])
		#ref_cluster_max <- rowMax(cnv[, idx_ref_cells])
		
		ref_cluster_min <- t(apply(cnv[, idx_ref_cells], 1, function(x) boxplot.stats(x)$stats[c(1,5)])) #lo whisker and high whisker
		ref_cluster_max <- ref_cluster_min[, 2]
		ref_cluster_min <- ref_cluster_min[, 1]
		
		temp_higher <- cnv
		idx_higher_tozero <- apply(temp_higher, 2, function(x) x < (ref_cluster_max + eps))
		temp_higher <- apply(temp_higher, 2, function(x) x - (ref_cluster_max + eps))
		temp_higher[idx_higher_tozero] <- 0
		
		temp_lower <- cnv
		idx_lower_tozero <- apply(temp_lower, 2, function(x) x > (ref_cluster_min - eps))
		temp_lower <- apply(temp_lower, 2, function(x) x - (ref_cluster_min - eps))
		temp_lower[idx_lower_tozero] <- 0
		
		temp_zero <- cnv
		idx_zero <-  apply(cnv, 2, function(x) ((ref_cluster_min - eps) < x) & (x < (ref_cluster_max + eps)))
		temp_zero[idx_zero] <- 0
		
		print(all(idx_higher_tozero | idx_lower_tozero | idx_zero))
		
		cnv <- temp_higher + temp_lower
	}
	
	return(list(cnv=cnv, ref_cluster=ref_cluster))
	
}
