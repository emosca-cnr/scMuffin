#' Adds a clustering to the clustering list
#' @param scMuffinList scMuffinList object
#' @param clusters named factor where the values are cells names and names are cluster labels
#' @param partition_id unique identifier for this partition
#' @description Add this partition to the set of partitions
#' @return scMuffinList with the additional partition
#' @importFrom stats setNames
#' @export

add_partitions <- function(scMuffinList=NULL, clusters=NULL, partition_id=NULL){
	
  
  clusters <- clusters[names(clusters) %in% colnames(scMuffinList$normalized)]
  
  if(length(clusters) < ncol(scMuffinList$normalized)){
    
    to_add <- which(!colnames(scMuffinList$normalized) %in% names(clusters))
    clusters <- c(clusters, setNames(rep(NA, length(to_add)), colnames(scMuffinList$normalized)[to_add]))
    clusters <- factor(clusters[match(colnames(scMuffinList$normalized), names(clusters))])
    clusters <- addNA(clusters)

  }
  
  if(length(scMuffinList$partitions)==0){
    
    cat("Creating 'scMuffinList$partitions'\n")
    scMuffinList$partitions <- setNames(data.frame(clusters, row.names = names(clusters)), partition_id)
    
  }else{
    
    idx_exist_part <- which(colnames(scMuffinList$partitions) == partition_id)
    
    if(length(idx_exist_part) > 0){ #overwrites
      
      cat("Overwriting an existing partition\n")
      
      scMuffinList$partitions[, idx_exist_part] <- clusters[match(rownames(scMuffinList$partitions), names(clusters))]
    
    }else{ #merge
      
      cat("Adding the partition\n")
      
      scMuffinList$partitions <- merge(scMuffinList$partitions, setNames(data.frame(clusters, row.names = names(clusters)), partition_id), by=0, sort=F, all=T)
      rownames(scMuffinList$partitions) <- scMuffinList$partitions$Row.names
      scMuffinList$partitions$Row.names <- NULL
      
    }
  }
  
	return(scMuffinList)
	
}