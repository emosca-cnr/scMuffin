#' Retrieve genes in detected CNV regions
#' 
#' @export
#' 
genes_in_detected_CNV_regions <- function(scMuffinList=NULL){
  
  ans <- CNV_PJ016$CNV$full$detected_cnv_regions
  for(i in 1:length(ans)){
    cat(names(ans)[i], "\n")
    df_i <- ans[[i]]
    ans[[i]] <- setNames(vector("list", nrow(df_i)), paste0(df_i$start, "____", df_i$stop))
    
    for(j in 1:nrow(df_i)){
      cat(j)
      reg_start <- which(rownames(CNV_PJ016$CNV$full$CNV)==df_i$start[j])
      reg_stop <- which(rownames(CNV_PJ016$CNV$full$CNV)==df_i$stop[j])
      
      CNV_j <- CNV_PJ016$CNV$full$CNV[reg_start:reg_stop, ]
      
      ans[[i]][[j]] <- regions_to_genes(CNV_j, CNV_PJ016$CNV$full$CNV_input)
      
      ans[[i]][[j]] <- unique(do.call(rbind, ans[[i]][[j]]))
      ans[[i]][[j]] <- ans[[i]][[j]][order(as.numeric(ans[[i]][[j]]$location), ans[[i]][[j]]$symbol), ]
      rownames(ans[[i]][[j]]) <- ans[[i]][[j]]$symbol
      
    }
  }
  
  return(ans)

}
