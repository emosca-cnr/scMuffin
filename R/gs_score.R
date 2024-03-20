#' Gene set scoring
#' @description Calculate gene set scores using the approach described in Tirosh2016
#' @param gene_set gene set
#' @param genes_by_cells real data
#' @param bins random data
#' @param nmark_min min number of markers
#' @param ncells_min min number of cells
#' @param k number of random gene set to be considered
#' @param kmin minimim number of k
#' @param verbose verbosity
#' @param null_model boolean, whether to use the permutations or not
#' @param na.rm whether to remove or not the missing values
#' @importFrom Matrix colMeans colSums
#' @return A data.frame with the following columns: 
#` \itemize{
#`  \item case: the gene set average observed in the cell; 
#`  \item case.N: number of genes in the gene set; 
#`  \item case.AV: number of genes available in the cell; 
#`  \item nmark_min: whether case.AV > nmark_min; 
#`  \item avg_control: average of the averag gene set score in control cells; 
#`  \item control.AV: number of available control cells;
#`  \item null_ok: whether a null could be defined for the cell; 
#`  \item avg_delta_score: average of case minu each control cell; 
#` }
#' @export
#' @references Tirosh2016 10.1126/science.aad0501


#gene_set_score <- function(gene_set_data, control_set_data, nmark_min=5, ncells_min=NULL){
gs_score <- function(gene_set=NULL, genes_by_cells=NULL, bins=NULL, nmark_min=5, ncells_min=5, k=100, kmin=50, verbose=FALSE, null_model=TRUE, na.rm=TRUE){
  
  if(sum(rownames(genes_by_cells) %in% gene_set) < nmark_min){
    
    message("not enough genes in gene set, please consider only gene set with enough number of genes.\n")
    res <- data.frame(case=rep(NA, ncol(genes_by_cells)), case.N=NA, case.AV=NA, nmark_min=F, avg_control=NA, control.AV=NA, null_ok=F, avg_delta_score=NA, stringsAsFactors = F, row.names = colnames(genes_by_cells))
    
  }else{
    
    cat("# gene set genes available:", sum(rownames(genes_by_cells) %in% gene_set), "/", length(gene_set), "\n")
    
    gene_set_data <- genes_by_cells[rownames(genes_by_cells) %in% gene_set, , drop=F]
    
    if(na.rm){
      gene_set_data[gene_set_data==0] <- NA #NA to avoid inclusion in mean calculation
    }
    
    if(null_model){
      control_set_data <- sc_create_null(genes_by_cells, bins = bins, gene_set = gene_set, k = k)
      if(verbose){
        message("nulls obtained: ", length(control_set_data), "\n")
      }
      for(i in 1:length(control_set_data)){
        if(na.rm){
          control_set_data[[i]][control_set_data[[i]]==0] <- NA
        }
      }
      
    }
    
    
    #new version
    ans <- list(gs=data.frame(case=colMeans(gene_set_data, na.rm = na.rm), case.N=nrow(gene_set_data), case.AV=colSums(!is.na(gene_set_data)), stringsAsFactors = F)) #real
    
    
    #attach to ans the genes-by-cells matrix of controls
    if(null_model){
      ans <- c(ans, lapply(control_set_data, function(x) data.frame(control=colMeans(x, na.rm = na.rm), control.N=nrow(x), control.AV=colSums(!is.na(x)), stringsAsFactors = F))) #controls
    }
    
    #case and control will have NaN for cells with all NA
    for(i in 1:length(ans)){
      ans[[i]][is.nan(ans[[i]][, 1]), 1] <- NA
    }
    
    ans <- lapply(ans, function(x) cbind(x, nmark_min=x[,3]>=nmark_min)) #add nmark_min flag
    
    if(null_model){
      
      #cells with enough genes
      vect_selected <- do.call(cbind, t(lapply(ans, function(x) x$nmark_min))) #the first column contains real data
      
      #null_ok <- rowSums(vect_selected[, -1])  > (ncol(vect_selected[, -1])) #old
      gene_set_ok_in_nulls <- rowSums(vect_selected[, -1], na.rm = na.rm) #number of permutations in which each cell satistfies nmark_min
      
      if(verbose){
        cat("gene_set_ok_in_nulls...\n")
        print(table(gene_set_ok_in_nulls[gene_set_ok_in_nulls>0]))
      }
      
      if(max(gene_set_ok_in_nulls) < kmin){
        cat("kmin", kmin, "; found, ", max(gene_set_ok_in_nulls), "\n")
        cat("Consider reducing kmin\n")
      }
      
      null_ok <- gene_set_ok_in_nulls >= kmin
      cat("# cells with valid permutations (at least nmark_min markers):", sum(null_ok), "\n")
      cat("# max k:", max(gene_set_ok_in_nulls), "\n")
      
      if(sum(vect_selected[, 1]) >= ncells_min){ #enough valid cells
        
        vect_selected[, 1] <- vect_selected[, 1] & null_ok #cells with enough valid permutations
        #vect_selected <- vect_selected & null_ok #cells with enough valid permutations
        
        if(sum(vect_selected[, 1]) >= ncells_min){
          
          #score matrix
          score <- do.call(cbind, t(lapply(ans[-1], function(x) ans[[1]]$case - x$control)))
          score[!vect_selected[, -1]] <- NA #differences without valid permutations (on the basis of nmark_min) are set to NA
          score <- rowMeans(score, na.rm = na.rm)
          score[is.nan(score)] <- NA
          
          #average control value for each cell
          avg_control <- do.call(cbind, t(lapply(ans[-1], function(x) x$control)))
          avg_control[!vect_selected[, -1]] <- NA #controls without valid permutations (on the basis of nmark_min) are set to NA
          avg_control <- rowMeans(avg_control, na.rm = na.rm)
          avg_control[is.nan(avg_control)] <- NA
          
          res <- data.frame(ans[[1]], avg_control=avg_control, control.AV=gene_set_ok_in_nulls, null_ok=null_ok, avg_delta_score=score, stringsAsFactors = F) #difference between real gs score (1), and each of the null values [2, k]
          
        }else{
          message("cannot create null models for at least ", ncells_min, " cells\n")
          res <- data.frame(ans[[1]], avg_control=NA, control.AV=NA, null_ok=F, avg_delta_score=NA, stringsAsFactors = F)
        }
        
        #cells without a sufficient number of genes
        res$avg_delta_score[!res$nmark_min | !res$null_ok] <- NA 

      }else{
        
        message("cannot calcolate gene set score for at least ", ncells_min, " cells\n")
        res <- data.frame(case=rep(NA, ncol(genes_by_cells)), case.N=NA, case.AV=NA, nmark_min=F, avg_control=NA, control.AV=NA, null_ok=F, avg_delta_score=NA, stringsAsFactors = F, row.names = colnames(genes_by_cells))
        
      }
      
    }else{ #without nullmodel
      
      res <- data.frame(ans[[1]], avg_control=NA, control.AV=NA, null_ok=NA, avg_delta_score=NA, stringsAsFactors = F)
      #cells without a sufficient number of genes
      res$case[!res$nmark_min] <- 0
    }
    
    #ans <- list(score_table=res, permutations=ans[-1])
    
  }
  
  return(res)
  
}
