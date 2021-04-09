#' merge_clusterings_enrichment
#' @export

merge_clusterings_enrichment <- function(clust_enrich){
	
	check_colnames <- lapply(clust_enrich, function(x) colnames(x$nes))
	for(i in 2:length(check_colnames)){
		if(!identical(check_colnames[[1]], check_colnames[[i]])){
			stop("not all colnames are identical\n")
		}
	}
	
	temp_nes <- clust_enrich[[1]]$nes
	rownames(temp_nes) <- paste0(names(clust_enrich)[1], "_", rownames(temp_nes))
	temp_fdrq <- clust_enrich[[1]]$fdrq
	rownames(temp_fdrq) <- paste0(names(clust_enrich)[1], "_", rownames(temp_fdrq))
	
	ans <- list(nes=temp_nes, fdrq=temp_fdrq)
	
	for(i in 2:length(clust_enrich)){
		
		temp_nes <- clust_enrich[[i]]$nes
		rownames(temp_nes) <- paste0(names(clust_enrich)[i], "_", rownames(temp_nes))
		
		temp_fdrq <- clust_enrich[[i]]$fdrq
		rownames(temp_fdrq) <- paste0(names(clust_enrich)[i], "_", rownames(temp_fdrq))
		
		ans$nes <- rbind(ans$nes, temp_nes[, match(colnames(ans$nes), colnames(temp_nes))])
		ans$fdrq <- rbind(ans$fdrq, temp_fdrq[, match(colnames(ans$fdrq), colnames(temp_fdrq))])
		
	}
	
	return(ans)
	
}