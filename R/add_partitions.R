#' Adds a clustering to the clustering list
#' @param scMuffinList clusterings object
#' @param clusters named factor with cluster label for each cell 
#' @param partition_id partition id
#' @description Adds a clustering to the parition list
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