#' prepare the gene set lists collected from various sources
#' @param gs_sources character vector, possible values: "SIG_CM_cancer", "SIG_CM_normal", "SIG_CancerSEA", "SIG_PNDB" and "msigdb"
#' @param custom_gsls named list of gene sets
#' @param CM_tissues character vector with strings to filter CellMArker db by tissue name
#' @param PNDB_tissues character vector with strings to filter Panglao db by tissue name
#' @param msigdb_hs_cat_subcat data.frame with species, category and subcategory: e.g. "Homo sapiens", "NA", "CP:KEGG"
#' @param genes_min minimum number of genes required in a gene set
#' @param genes_max maximum number of genes required in a gene set
#' @param genes optional character vectors with the universe of all the genes under analysis
#' @importFrom utils data
#' @export

prepare_gsls <- function(gs_sources=NULL, custom_gsls=NULL, CM_tissues=NULL, PNDB_tissues=NULL, msigdb_hs_cat_subcat=NULL, genes_min=5, genes_max=500, genes=NULL){
	
	SIG_CM_normal <- SIG_CM_cancer <- SIG_CancerSEA <- SIG_PNDB <- NULL #to please the check
	gsls <- custom_gsls
	
	if(!is.null(gs_sources)){
		
		if(!is.null(CM_tissues)){
			CM_tissues <- paste("CM__", CM_tissues, "__", sep = "")
		}
		
		if(!is.null(PNDB_tissues)){
			PNDB_tissues <- paste("PN__", PNDB_tissues, "__", sep = "")
		}
		
		if(!is.null(custom_gsls)){
			gsls <- custom_gsls
		}else{
			gsls <- list()
		}
		
		if("SIG_CM_cancer" %in% gs_sources){
			data("SIG_CM_cancer", envir=environment())
			temp <- list()
			for(i in 1:length(CM_tissues)){
				temp <- c(temp, (SIG_CM_cancer[grepl(paste0("^", CM_tissues[i]), names(SIG_CM_cancer))]))
			}
			gsls$SIG_CM_cancer <- temp
		}
		if("SIG_CM_normal" %in% gs_sources){
			data("SIG_CM_normal", envir=environment())
			temp <- list()
			for(i in 1:length(CM_tissues)){
				temp <- c(temp, (SIG_CM_normal[grepl(paste0("^", CM_tissues[i]), names(SIG_CM_normal))]))
			}
			gsls$SIG_CM_normal <- temp
		}
		if("SIG_CancerSEA" %in% gs_sources){
			data("SIG_CancerSEA", envir=environment())
			gsls$SIG_CancerSEA <- SIG_CancerSEA
		}
		if("SIG_PNDB" %in% gs_sources){
			data("SIG_PNDB", envir=environment())
			temp <- list()
			for(i in 1:length(PNDB_tissues)){
				temp <- c(temp, (SIG_PNDB[grepl(paste0("^", PNDB_tissues[i]), names(SIG_PNDB))]))
			}
			gsls$SIG_PNDB <- temp
		}
		
		if("msigdb" %in% gs_sources & !is.null(msigdb_hs_cat_subcat)){
			msigdb_sigs <- apply(msigdb_hs_cat_subcat, 1, function(x) paste(x, collapse = "_"))
			for(i in 1:nrow(msigdb_hs_cat_subcat)){
				gsls[[msigdb_sigs[i]]] <- get_msigdb_geneset(species = msigdb_hs_cat_subcat[i, 1], category = msigdb_hs_cat_subcat[i, 2], subcategory = msigdb_hs_cat_subcat[i, 3])$path_list
			}
		}
		
	}
	
	if(!is.list(gsls)){
	  message("Empty gene set list\n")
	}
	
	### filter gsls
	if(!is.null(genes)){
		gsls <- lapply(gsls, function(x) lapply(x, function(y) unique(y[y %in% genes])))
	}
	for(i in 1:length(gsls)){
		gs_length <- unlist(lapply(gsls[[i]], length))
		gsls[[i]] <- gsls[[i]][gs_length >= genes_min & gs_length <= genes_max]
	}
	
	#gsls_length <- unlist(lapply(gsls[[i]], length))
	
	
	return(gsls)
	
}