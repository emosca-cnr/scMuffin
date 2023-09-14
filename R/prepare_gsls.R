#' Prepare the gene set lists collected from various sources
#' @description Prepare the gene set lists collected from various sources
#' @param gs_sources character vector, possible values: "CM_cancer", "CM_normal", "CancerSEA", "PNDB" and "msigdb"
#' @param custom_gsls named list of gene sets
#' @param CM_tissues character vector with strings to filter CellMArker db by tissue name
#' @param PNDB_tissues character vector with strings to filter Panglao db by tissue name
#' @param msigdb_hs_cat_subcat data.frame with species, category and subcategory: e.g. "Homo sapiens", "NA", "CP:KEGG"
#' @param genes_min minimum number of genes required in a gene set
#' @param genes_max maximum number of genes required in a gene set
#' @param scMuffinList scMuffinList object
#' @param id_type type of id: Symbol or EntrezID
#' @param genes all genes to be considered
#' @importFrom utils data
#' @return A lists of gene set lists
#' @export

prepare_gsls <- function(gs_sources=NULL, custom_gsls=NULL, CM_tissues=NULL, PNDB_tissues=NULL, msigdb_hs_cat_subcat=NULL, genes_min=5, genes_max=500, id_type=c("Symbol", "EntrezID"), scMuffinList=NULL, genes=NULL){
  
  #SIG_CM_normal <- SIG_CM_cancer <- SIG_CancerSEA <- SIG_PNDB <- NULL #to please the check
  gsls_EntrezID <- gsls_Symbol <- NULL #to please the check
  gsls <- custom_gsls
  
  if(is.null(genes)){
    if(length(scMuffinList$normalized)==0){
      stop("scMuffinList does not contain a normalized genes_by_cells expression matrix\n")
    }
  }
  
  id_type <- match.arg(id_type, c("Symbol", "EntrezID"))
  
  if(id_type == "Symbol"){
    
    data("gsls_Symbol", envir=environment())
    gsls_source <- gsls_Symbol
    
  }else{
    
    data("gsls_EntrezID", envir=environment())
    gsls_source <- gsls_EntrezID
    
  }
  
  if(!is.null(gs_sources)){
    
    if(!is.null(CM_tissues)){
      CM_tissues <- paste("CM__", CM_tissues, "__", sep = "") 
    }
    
    if(!is.null(PNDB_tissues)){
      PNDB_tissues <- paste("PN__", PNDB_tissues, "__", sep = "")
    }
    
    if(!is.null(custom_gsls)){
      if(!is.list(custom_gsls)){
        stop("custom_gsls must be a list\n.")
      }
      gsls <- custom_gsls
    }else{
      gsls <- list()
    }
    
    if("CM_cancer" %in% gs_sources){
      #data("SIG_CM_cancer", envir=environment())
      temp <- list()
      for(i in 1:length(CM_tissues)){
        #temp <- c(temp, (SIG_CM_cancer[grepl(paste0("^", CM_tissues[i]), names(SIG_CM_cancer))]))
        temp <- c(temp, (gsls_source$CM_cancer[grepl(paste0("^", CM_tissues[i]), names(gsls_source$CM_cancer))]))
      }
      if(length(temp)==0){
        cat("None of the given CM_tissues was found among CM_cancer gene sets.\n")
      }
      gsls$CM_cancer <- temp
    }
    
    if("CM_normal" %in% gs_sources){
      #data("SIG_CM_normal", envir=environment())
      temp <- list()
      for(i in 1:length(CM_tissues)){
        #temp <- c(temp, (SIG_CM_normal[grepl(paste0("^", CM_tissues[i]), names(SIG_CM_normal))]))
        temp <- c(temp, (gsls_source$CM_normal[grepl(paste0("^", CM_tissues[i]), names(gsls_source$CM_normal))]))
      }
      if(length(temp)==0){
        cat("None of the given CM_tissues was found among CM_normal gene sets.\n")
      }
      gsls$CM_normal <- temp
    }
    
    if("CancerSEA" %in% gs_sources){
       gsls$CancerSEA <- gsls_source$CancerSEA
    }
    
    if("PNDB" %in% gs_sources){
      #data("SIG_PNDB", envir=environment())
      temp <- list()
      for(i in 1:length(PNDB_tissues)){
        #temp <- c(temp, (SIG_PNDB[grepl(paste0("^", PNDB_tissues[i]), names(SIG_PNDB))]))
        temp <- c(temp, (gsls_source$PNDB[grepl(paste0("^", PNDB_tissues[i]), names(gsls_source$PNDB))]))
      }
      if(length(temp)==0){
        cat("None of the given CM_tissues was found among PNDB gene sets.\n")
      }
      
      gsls$PNDB <- temp
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
  if(is.null(genes)){
    genes <- rownames(scMuffinList$normalized)
    idx_zero <- which(rowSums(scMuffinList$normalized)==0)
    if(length(idx_zero)>0){
      cat("Found genes with all-zero values. These genes may cause issues. Trying to proceed removing these genes from the gene sets.\n")
      genes <- genes[-idx_zero]
    }
  }
  gsls <- lapply(gsls, function(x) lapply(x, function(y) unique(y[y %in% genes])))
  
  
  cat("Current gene set list size.\n")
  print(lapply(gsls, lengths))
  
  cat("Filtering: [", genes_min, ",", genes_max, "].\n")
  for(i in 1:length(gsls)){
    gs_length <- unlist(lapply(gsls[[i]], length))
    gsls[[i]] <- gsls[[i]][gs_length >= genes_min & gs_length <= genes_max]
  }

  cat("Final gene set list size.\n")
  print(lapply(gsls, lengths))
  
  return(gsls)
  
}