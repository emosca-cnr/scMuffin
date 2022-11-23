#' Add a feature to scMuffinList
#' @param scMuffinList scMuffinList object
#' @param name feature name
#' @param summary data.frame with cell identifiers as row.names
#' @param full optional information related to the feature, like cell-level statistics
#' @return scMuffinList with the additional feature
#' @export

add_features <- function(scMuffinList=NULL, name=NULL, summary=NULL, full=NULL){
	
  scMuffinList[[name]] <- list(summary=summary, full=full)

	return(scMuffinList)
	
}