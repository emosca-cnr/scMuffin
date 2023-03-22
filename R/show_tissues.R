#' Show the tissues available for \code{\link{prepare_gsls}}
#' @param gsc_source the source database
#' @description Show the tissues available by means of \code{\link{prepare_gsls}}
#' @return List of gsc_source with tissue vectors.
#' @export

show_tissues <- function(gsc_source=c("PNDB", "CM_normal", "CM_cancer")){
	
	ans <- list()
	gsls_EntrezID <- gsls_Symbol <- NULL #to please the check
	data("gsls_Symbol", envir=environment())
	
	if("PNDB" %in% gsc_source){
		ans$PNDB <- unique(gsub("^PN__(.+)__.+$", "\\1", names(gsls_Symbol$PNDB)))
	}
	
	if("CM_normal" %in% gsc_source){
	  ans$CM_normal <- unique(gsub("^CM__(.+)__Normal__.+$", "\\1", names(gsls_Symbol$CM_normal)))
	}
	
	if("CM_cancer" %in% gsc_source){
	  ans$CM_cancer <- unique(gsub("^CM__(.+)__.+__.+$", "\\1", names(gsls_Symbol$CM_cancer)))
	}
	
	return(ans)
	
}


