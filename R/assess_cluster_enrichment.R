#' assess_cluster_enrichment
#' @param features feature list
#' @param clusterings clustering list
#' @param meta_clusters boolean, TRUE if the assessment involves meta clusters
#' @param write_output boolean, TRUE to write the output in the output dir
#' @param out_dir output dir
#' @export
#' @importFrom utils write.table
assess_cluster_enrichment <- function(features, clusterings, meta_clusters=FALSE, write_output=TRUE, out_dir="./"){
  
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive=TRUE)
  }
  
  X_fact <- as.matrix(features$df[, features$type=="factor"])
  X_num <- as.matrix(features$df[, !features$type=="factor"])
  X_num[is.na(X_num)] <- 0
  
  
  if(meta_clusters){
    
    ##### NUMERCIC -> GSEA
    if(ncol(X_num)>0){
      
      cluster_gsea_res_nes <- cluster_gsea(X_num, setNames(clusterings$meta_cl, clusterings$cell_id))
      
      if(write_output){
        for(j in 1:nrow(cluster_gsea_res_nes$fdrq)){
          out_table <- data.frame(feature=colnames(cluster_gsea_res_nes$nes), nes=cluster_gsea_res_nes$nes[j, ], fdrq=cluster_gsea_res_nes$fdrq[j, ], stringsAsFactors = F)
          write.table(out_table, row.names = F, sep="\t", file = paste0(out_dir, "/cluster_enrichment_", rownames(cluster_gsea_res_nes$fdrq)[j], ".txt"))
        }
      }
    }
    
    ##### FACTOR -> Fisher exact test
    if(ncol(X_fact)>0){
      
    }
    
  }else{
    
    ##### NUMERCIC -> GSEA
    if(ncol(X_num)>0){
      
      cluster_gsea_res_nes <- vector("list", ncol(clusterings))
      names(cluster_gsea_res_nes) <- colnames(clusterings)
      cluster_gsea_res_fdr <- cluster_gsea_res_nes
      
      for(i in 1:ncol(clusterings)){
        cat(i)
        cluster_gsea_res_nes[[i]] <- cluster_gsea(X_num, setNames(clusterings[, i], rownames(clusterings)))
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
    
    
    ##### FACTOR -> Fisher exact test
    
    
    
  }
  
  
  return(cluster_gsea_res_nes)
  
}