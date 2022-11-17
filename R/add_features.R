#' Adds a feature to the feature list
#' @param features_obj features object
#' @param to_add data.frame of 1 or more cell values with cell id as row names
#' @description Adds a feature to the feature list
#' @return The scMuffinList with the the additional feature
#' @export

add_features <- function(scMuffinList=NULL, id=NULL, summary=NULL, full=NULL){
	
  
  scMuffinList[[id]] <- list(summary=summary, full=full)
  

	return(scMuffinList)
	
}