#' assess_cluster_enrichment
#' @export
#' @importFrom utils write.table
assess_cluster_enrichment <- function(features, clusterings, meta_clusters=FALSE, write_output=TRUE, out_dir="./"){
	
	if(!dir.exists(out_dir)){
		dir.create(out_dir, recursive=TRUE)
	}
	if(meta_clusters){
		
		cluster_gsea_res_nes <- cluster_gsea(features, setNames(clusterings$meta_cl, clusterings$cell_id))
		
		if(write_output){
			for(j in 1:nrow(cluster_gsea_res_nes$fdrq)){
				out_table <- data.frame(feature=colnames(cluster_gsea_res_nes$nes), nes=cluster_gsea_res_nes$nes[j, ], fdrq=cluster_gsea_res_nes$fdrq[j, ], stringsAsFactors = F)
				write.table(out_table, row.names = F, sep="\t", file = paste0(out_dir, "/cluster_enrichment_", rownames(cluster_gsea_res_nes$fdrq)[j], ".txt"))
			}
		}
		
	}else{
		
		cluster_gsea_res_nes <- vector("list", ncol(clusterings))
		names(cluster_gsea_res_nes) <- colnames(clusterings)
		cluster_gsea_res_fdr <- cluster_gsea_res_nes
		
		for(i in 1:ncol(clusterings)){
			cat(i)
			cluster_gsea_res_nes[[i]] <- cluster_gsea(features, setNames(clusterings[, i], rownames(clusterings)))
		}
		
		if(write_output){
			for(i in 1:length(cluster_gsea_res_nes)){
				for(j in 1:nrow(cluster_gsea_res_nes[[i]]$fdrq)){
					
					out_table <- data.frame(feature=colnames(cluster_gsea_res_nes[[i]]$nes), nes=cluster_gsea_res_nes[[i]]$nes[j, ], fdrq=cluster_gsea_res_nes[[i]]$fdrq[j, ], stringsAsFactors = F)
					write.table(out_table, row.names = T, sep="\t", file = paste0(out_dir, "/ce_", names(cluster_gsea_res_nes)[i], "_", rownames(cluster_gsea_res_nes[[i]]$fdrq)[j], ".txt"))
				}
			}
			
		}
		
	}
	
	
	return(cluster_gsea_res_nes)
	
}