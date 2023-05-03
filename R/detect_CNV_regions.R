#' Detect CNV regions
#' @param scMuffinList scMuffinList object
#' @param z.score whether to use clustere median z-scores of CNV signal (TRUE) or cluster median CNV signal (FALSE).
#' @param eps the value that defines a CNV. If NUll the stard deviation of the dataset is used.
#' @description Plot an heatmap of the CNV.
#' @details CNV Profile of every cluster
#' @import ComplexHeatmap grDevices grid
#' @export

detect_CNV_regions <- function(scMuffinList = NULL, z.score=FALSE, eps=NULL){
  
  if(length(scMuffinList$partitions[, "CNV"])==0){
    stop("CNV regions detection requires scMuffinList$partitions[, 'CNV']\n")
  }
  
  cnv_clustering <- scMuffinList$partitions[, "CNV"]
  names(cnv_clustering) <- rownames(scMuffinList$partitions)
  cnv <- scMuffinList$CNV$full$CNV
  
  cnv <- cnv[, match(names(cnv_clustering), colnames(cnv))]
  
  row_chr <- row_chr_vector <- gsub("(chr[^_]+)_.+", "\\1", rownames(cnv))
  row_chr <- factor(row_chr, levels = unique(row_chr))
  row_chr <- split(row_chr, row_chr)
  ngenes_chr <- unlist(lapply(row_chr, length))
  
  if(z.score){
    avg_cl_cnv <- t(apply(t(scale(t(cnv))), 1, function(i_row) tapply(i_row, cnv_clustering, median)))
  }else{
    avg_cl_cnv <- t(apply(cnv, 1, function(i_row) tapply(i_row, cnv_clustering, median)))
  }
  if(is.null(eps)){
    eps <- sd(as.numeric(avg_cl_cnv))
  }
  cat("Using eps =", eps, "\n")
  
  cnv_detected <- sign(abs(avg_cl_cnv)>eps)
  cnv_detected <- split.data.frame(cnv_detected, row_chr_vector)
  cnv_begin_end <- lapply(cnv_detected, function(x) x[-1, , drop=FALSE] - x[-nrow(x), ,drop=FALSE])
  cnv_begin_end <- lapply(cnv_begin_end, function(x) x[rowSums(abs(x))>0, , drop=F])
  
  #adjust the beginning of every chromosome
  for(i in 1:length(cnv_begin_end)){
    
    ### add zero's at the begin and at the end
    cnv_begin_end[[i]] <- rbind(
      as.data.frame(matrix(0, nrow=1, ncol = ncol(cnv_begin_end[[i]]), dimnames = list(rownames(cnv_detected[[i]])[1], colnames(cnv_begin_end[[i]])))),
      cnv_begin_end[[i]])
    
    for(j in 1:ncol(cnv_begin_end[[i]])){
      
      #start stop values
      temp <- cnv_begin_end[[i]][cnv_begin_end[[i]][, j]!=0, j]
      
      if(length(temp)==0 & all(cnv_detected[[i]][, j]==1)){
        cnv_begin_end[[i]][1, j] <- 1
        cnv_begin_end[[i]][nrow(cnv_begin_end[[i]]), j] <- -1
      }
      
      if(length(temp)>0){
        
        
        firstval <- temp[1]
        lasttval <- temp[length(temp)]
        
        if(firstval==-1){
          cnv_begin_end[[i]][1, j] <- 1
        }
        if(lasttval==1){
          cnv_begin_end[[i]][nrow(cnv_begin_end[[i]]), j] <- -1
        }
      }
    }
  }
  
  cnv_regions <- setNames(vector("list", length = length(cnv_begin_end)), names(cnv_begin_end))
  for(i in 1:length(cnv_begin_end)){
    for(j in 1:ncol(cnv_begin_end[[i]])){
      temp <- cnv_begin_end[[i]][cnv_begin_end[[i]][, j] !=0, j, drop=F]
      temp_l <- nrow(temp)
      
      if(temp_l>0){
        temp <- data.frame(region=rownames(temp), pattern=temp[,1])
        
        #if it ends with two -1 or is just a 1 -> add two lines: 1 and -1
        if(temp$pattern[temp_l] == -1){
          if(temp_l>1){
            if(temp$pattern[temp_l-1]==-1){
              temp$pattern[temp_l] <- 1
              temp <- rbind(temp, data.frame(region=temp$region[temp_l], pattern=-1))
            }
          }else{
            temp$pattern[temp_l] <- 1
            temp <- rbind(temp, data.frame(region=temp$region[temp_l], pattern=-1))
          }
        }
        
        temp <- data.frame(chr=gsub("^(chr[^_]+)__.+$", "\\1", temp$region[temp$pattern==1]), start=temp$region[temp$pattern==1], start.loc=as.numeric(gsub("^.+_(\\d+)__.+$", "\\1", temp$region[temp$pattern==1])), stop=temp$region[temp$pattern==-1], stop.loc=as.numeric(gsub("^.+_(\\d+)$", "\\1", temp$region[temp$pattern==-1])), cluster=colnames(cnv_begin_end[[i]])[j], stringsAsFactors = F)
        temp$length <- temp$stop.loc - temp$start.loc
        
        if(length(cnv_regions[[i]])==0){
          cnv_regions[[i]] <- temp
        }else{
          cnv_regions[[i]] <- rbind(cnv_regions[[i]], temp)
        }
      }
    }
  }
  
  cnv_regions <- cnv_regions[unlist(lapply(cnv_regions, length))>0]
  
  #concatenete adjcent regions
  ans <- cnv_regions
  for(i in 1:length(ans)){
    ans[[i]] <- ans[[i]][1, , drop=F]
    if(nrow(cnv_regions[[i]])>1){
      ans_i <- 1
      for(j in 2:nrow(cnv_regions[[i]])){
        #if the j starts before j-1 finish replace the stop of j-1 with that of j
        if(cnv_regions[[i]]$start.loc[j] < cnv_regions[[i]]$stop.loc[j-1] & cnv_regions[[i]]$cluster[j] == cnv_regions[[i]]$cluster[j-1]){
          ans[[i]]$stop.loc[ans_i] <- cnv_regions[[i]]$stop.loc[j]
          ans[[i]]$stop[ans_i] <- cnv_regions[[i]]$stop[j]
          ans[[i]]$length[ans_i] <- ans[[i]]$stop.loc[ans_i] - ans[[i]]$start.loc[ans_i]
        }else{
          ans_i <- ans_i + 1
          ans[[i]] <- rbind(ans[[i]], cnv_regions[[i]][j, , drop=F])
        }
      }
    }
  }
  
  
  return(ans)
  
}
