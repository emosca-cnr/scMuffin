#' Map lanscent ppi to symbols
#' @import org.Hs.eg.db
#' @author Ettore Mosca
map_landscent_ppi <- function(){
	
	ppi <- LandSCENT::net17Jan16.m
	eg2sym <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
	eg2sym <- eg2sym[eg2sym$gene_id %in% rownames(ppi), ]
	eg2sym <- merge(eg2sym, data.frame(gene_id=rownames(ppi), stringsAsFactors = F), by="gene_id", sort=F)
	
	if(length(unique(eg2sym$symbol)) == length(unique(eg2sym$gene_id))){
	
		rownames(ppi) <- colnames(ppi) <- eg2sym$symbol[match(rownames(ppi), eg2sym$gene_id)]
	
	}else{
		warning("mapping produced duplicates, returning the initial matrix\n")
	}

	return(ppi)

}