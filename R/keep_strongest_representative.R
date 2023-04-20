#' Keep the feature with the average highest values
#' @param genes_by_cells genes-by-cells expression matrix
#' @param row.names row.names that have to be assigned
#' @export
#' 
keep_strongest_representative <- function(genes_by_cells=NULL, row.names=NULL){
	
	rn_avg <- data.frame(rn=row.names, avg=rowMeans(genes_by_cells), genes_by_cells, stringsAsFactors = F)
	rn_avg <- rn_avg[order(-rn_avg$avg), ]
	idx_dupl <- which(duplicated(rn_avg$rn))
	rn_avg <- rn_avg[-idx_dupl, ]
	rownames(rn_avg) <- rn_avg$rn
	rn_avg$avg <- rn_avg$rn <- NULL
	rn_avg <- rn_avg[order(rownames(rn_avg)), ]
	
	return(rn_avg)
	
}