#' Returns msigdb geneset in a format compatible with calculate_signatures()
#' @param species Species name
#' @param category MSigDB category name, see msigdbr_collections()
#' @param subcategory MsigDB subcategory name, see msigdbr_collections()
#' @param type Gene name of interest, can be gene_symbol or entrez_gene
#' @description Returns msigdb geneset in a format compatible with calculate_signatures()
#' @import msigdbr
#' @return List with msigdb_output, the output of msigdbr::msigdbr, and path_list, the list of gene sets
#' @export

get_msigdb_geneset <- function(species, category=NULL, subcategory=NULL, type="gene_symbol") {
	if(is.na(category)){
		category <- NULL
	}
	if(is.na(subcategory)){
		subcategory <- NULL
	}
  msig <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
  if(type=="gene_symbol") {
    msig_list <- split(x=msig$gene_symbol, f = msig$gs_name)
  } else if (type == "entrez_gene") {
    msig_list <- split(x=msig$entrez_gene, f = msig$gs_name)
  }
  output <- list(msigdb_output = msig, path_list = msig_list)
  return(output)
}