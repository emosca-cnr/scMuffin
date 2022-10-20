#' Show the tissues available by means of prepare_gsl
#' @param gsc_source the source database
#' @description Show the tissues available by means of prepare_gsl
#' 
#' @export

show_tissues <- function(gsc_source=c("SIG_PNDB", "SIG_CM_normal", "SIG_CM_cancer")){
	
	ans <- list()
	SIG_CM_normal <- SIG_CM_cancer <- SIG_CancerSEA <- SIG_PNDB <- NULL #to please the check
	
	if("SIG_PNDB" %in% gsc_source){
	  data("SIG_PNDB", envir=environment())
		ans$SIG_PNDB <- unique(gsub("^PN__(.+)__.+$", "\\1", names(SIG_PNDB)))
	}
	
	if("SIG_CM_normal" %in% gsc_source){
	  data("SIG_CM_normal", envir=environment())
	  ans$SIG_CM_normal <- unique(gsub("^CM__(.+)__Normal__.+$", "\\1", names(SIG_CM_normal)))
	}
	
	if("SIG_CM_cancer" %in% gsc_source){
	  data("SIG_CM_cancer", envir=environment())
	  ans$SIG_CM_cancer <- unique(gsub("^CM__(.+)__.+__.+$", "\\1", names(SIG_CM_cancer)))
	}
	
	return(ans)
	
}


