#' cluster_stats
#' @param feat_obj feature object
#' @param partitions partitions object
#' @param which_col which column of the clustering object should be used. if NULL, all columns will be considered
#' @param avg function to use for average
#' @param var function to use for variability
#' @param na.rm whether to remove na or not
#' @importFrom stats sd
#' @export

cluster_stats <- function(feat_obj=NULL, partitions=NULL, which_col=NULL, avg=mean, var=sd, na.rm=TRUE){
	
	partitions <- partitions[match(rownames(feat_obj$df), rownames(partitions)), , drop=F]
		
	if(is.null(which_col)){
		ans_avg <- ans_var <- vector("list", ncol(partitions))
		names(ans_avg) <- names(ans_var) <- colnames(partitions)
		for(i in 1:length(ans_avg)){
			ans_avg[[i]] <- t(apply(feat_obj$df[, feat_obj$type == "numeric"], 2, function(x) tapply(x, partitions[, i], avg, na.rm=na.rm)))
			ans_avg[[i]][is.nan(ans_avg[[i]]) | is.na(ans_avg[[i]])] <- 0
			
			ans_var[[i]] <- t(apply(feat_obj$df[, feat_obj$type == "numeric"], 2, function(x) tapply(x, partitions[, i], var, na.rm=na.rm)))
			ans_var[[i]][is.nan(ans_var[[i]]) | is.na(ans_var[[i]])] <- 0
		}
	}else{
		ans_avg <- t(apply(feat_obj$df[, feat_obj$type == "numeric"], 1, function(x) tapply(x, partitions[, which_col], avg, na.rm)))
		ans_var <- t(apply(feat_obj$df[, feat_obj$type == "numeric"], 1, function(x) tapply(x, partitions[, which_col], var, na.rm)))
	}
	
	return(list(avg=ans_avg, var=ans_var))
	
}