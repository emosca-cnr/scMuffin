#' Adds a clustering to the clustering list
#' @param scMuffinList scMuffinList object
#' @param clusters named factor where the values are cells names and names are cluster labels
#' @param partition_id unique identifier for this partition
#' @description Add this partition to the set of partitions
#' @return scMuffinList with the additional partition
#' @export

add_partitions <- function(scMuffinList=NULL, clusters=NULL, partition_id=NULL){
	
  if(length(scMuffinList$partitions)==0){
    scMuffinList$partitions <- setNames(data.frame(clusters, row.names = names(clusters)), partition_id)
  }else{
    scMuffinList$partitions <- merge(scMuffinList$partitions, setNames(data.frame(clusters, row.names = names(clusters)), partition_id), by=0, sort=F, all=T)
    rownames(scMuffinList$partitions) <- scMuffinList$partitions$Row.names
    scMuffinList$partitions$Row.names <- NULL
  }
  
	return(scMuffinList)
	
}