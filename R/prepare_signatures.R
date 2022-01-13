#' prepare_signatures collected from various sources
#' @param signature_sources character vector, possible values: "SIG_CM_cancer", "SIG_CM_normal", "SIG_CSEA", "SIG_PNDB" and "msigdb"
#' @param custom_signatures named list of gene sets
#' @param CM_tissues character vector with strings to filter CellMArker db by tissue name
#' @param PNDB_tissues character vector with strings to filter Panglao db by tissue name
#' @param msigdb_hs_cat_subcat data.frame with species, category and subcategory: e.g. "Homo sapiens", "NA", "CP:KEGG"
#' @param genes_min minimum number of genes required in a gene set
#' @param genes_max maximum number of genes required in a gene set
#' @param genes optional character vectors with the universe of all the genes under analysis
#' @importFrom utils data
#' @export

prepare_signatures <- function(signature_sources=NULL, custom_signatures=NULL, CM_tissues=NULL, PNDB_tissues=NULL, msigdb_hs_cat_subcat=NULL, genes_min=5, genes_max=500, genes=NULL){
	
	SIG_CM_normal <- SIG_CM_cancer <- SIG_CSEA <- SIG_PNDB <- NULL #to please the check
	signatures <- custom_signatures
	
	if(!is.null(signature_sources)){
		
		if(!is.null(CM_tissues)){
			CM_tissues <- paste("CM__", CM_tissues, "__", sep = "")
		}
		
		if(!is.null(PNDB_tissues)){
			PNDB_tissues <- paste("PN__", PNDB_tissues, "__", sep = "")
		}
		
		if(!is.null(custom_signatures)){
			signatures <- custom_signatures
		}else{
			signatures <- list()
		}
		
		if("SIG_CM_cancer" %in% signature_sources){
			data("SIG_CM_cancer", envir=environment())
			temp <- list()
			for(i in 1:length(CM_tissues)){
				temp <- c(temp, (SIG_CM_cancer[grepl(paste0("^", CM_tissues[i]), names(SIG_CM_cancer))]))
			}
			signatures$SIG_CM_cancer <- temp
		}
		if("SIG_CM_normal" %in% signature_sources){
			data("SIG_CM_normal", envir=environment())
			temp <- list()
			for(i in 1:length(CM_tissues)){
				temp <- c(temp, (SIG_CM_normal[grepl(paste0("^", CM_tissues[i]), names(SIG_CM_normal))]))
			}
			signatures$SIG_CM_normal <- temp
		}
		if("SIG_CSEA" %in% signature_sources){
			data("SIG_CSEA", envir=environment())
			signatures$SIG_CSEA <- SIG_CSEA
		}
		if("SIG_PNDB" %in% signature_sources){
			data("SIG_PNDB", envir=environment())
			temp <- list()
			for(i in 1:length(PNDB_tissues)){
				temp <- c(temp, (SIG_PNDB[grepl(paste0("^", PNDB_tissues[i]), names(SIG_PNDB))]))
			}
			signatures$SIG_PNDB <- temp
		}
		
		if("msigdb" %in% signature_sources & !is.null(msigdb_hs_cat_subcat)){
			msigdb_sigs <- apply(msigdb_hs_cat_subcat, 1, function(x) paste(x, collapse = "_"))
			for(i in 1:nrow(msigdb_hs_cat_subcat)){
				signatures[[msigdb_sigs[i]]] <- get_msigdb_geneset(species = msigdb_hs_cat_subcat[i, 1], category = msigdb_hs_cat_subcat[i, 2], subcategory = msigdb_hs_cat_subcat[i, 3])$path_list
			}
		}
		
	}
	### filter signatures
	if(!is.null(genes)){
		signatures <- lapply(signatures, function(x) lapply(x, function(y) unique(y[y %in% genes])))
	}
	for(i in 1:length(signatures)){
		signatures_length <- unlist(lapply(signatures[[i]], length))
		signatures[[i]] <- signatures[[i]][signatures_length >= genes_min & signatures_length <= genes_max]
	}
	
	signatures_length <- unlist(lapply(signatures[[i]], length))
	
	
	return(signatures)
	
}