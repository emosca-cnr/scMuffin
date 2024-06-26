#' Keep the feature with the highest average across cells
#' @param genes_by_cells genes-by-cells expression matrix
#' @param row.names row.names that have to be assigned
#' @export
#' 
keep_strongest_representative <- function(genes_by_cells=NULL, row.names=NULL){
  
  #calculate rowMeans
  rn_avg <- data.frame(rn=row.names, avg=rowMeans(genes_by_cells), genes_by_cells, stringsAsFactors = F)
  
  #sort by decreasing rowMeans
  rn_avg <- rn_avg[order(-rn_avg$avg), ]
  
  #find duplicates among novel row.names
  idx_dupl <- which(duplicated(rn_avg$rn))
  if(length(idx_dupl)>0){
    rn_avg <- rn_avg[-idx_dupl, ]
    
    #if there are NA, remove the row
    idx_NA <- which(is.na(rn_avg$rn))
    if(length(idx_NA)>0){
      rn_avg <- rn_avg[-idx_NA, ]
    }
    
    rownames(rn_avg) <- rn_avg$rn
    rn_avg$avg <- rn_avg$rn <- NULL
    rn_avg <- rn_avg[order(rownames(rn_avg)), ]
  }else{
    cat("No duplicated features.\n")
    rn_avg <- NULL
  }
  
  return(rn_avg)
  
}