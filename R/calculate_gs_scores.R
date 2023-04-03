#' Calculate gene set scores
#' @description Calculate gene set scores using the approach described in Tirosh2016
#' @param scMuffinList scMuffinList object
#' @param gs_list list of gene sets
#' @param mc.cores number of cores
#' @param nbins number of bins to split the distribution of average gene expression
#' @param nmark_min number of minimum markers that are required for the succesful calculation of a gene set score
#' @param ncells_min number of minimum cells in which a gene set has to be succesfully calculated
#' @param k number of permutations
#' @param kmin minimum number of permutations; due to missing values it is hard to ensure that a gene set score can be compared to k permutations in every cell
#' @param score_type type of score. if "relative", than the score is the difference between the observed gene set average expression and that of a k permutations; if "mean" the score is equal to the observed gene set average expression
#' @param null_model TRUE if permutations have to be used. Required for score_type="relative"
#' @param verbose verbosity
#' @param na.rm whether to use NA or not
#' @param overwrite whether to update or not gene_set_scoring and gene_set_scoring_full elements of scMuffinList. 
#' @return scMuffinList with element gene_set_scoring, a list that contains summary and full. The element summary contains a cells-by-gene sets data.frame. The element "full" contains a data.frame for each gene set. See [gs_score()] for further details.
#' @references Tirosh2016 10.1126/science.aad0501
#' 
#' @export
#' @import parallel

calculate_gs_scores <- function(scMuffinList=NULL, gs_list=NULL, mc.cores=2, nbins=25, nmark_min = 5, ncells_min = 10, k=100, kmin=50, score_type=c("relative", "mean"), null_model=TRUE, verbose=TRUE, na.rm=TRUE, overwrite=FALSE){
  
  
  cat("nbins:", nbins, "\n")
  cat("nmark_min:", nmark_min, "\n")
  cat("ncells_min:", ncells_min, "\n")
  cat("k:", k, "\n")
  cat("kmin:", kmin, "\n")
  cat("score_type:", score_type, "\n")
  cat("null_model:", null_model, "\n")
  cat("verbose:", verbose, "\n")
  cat("na.rm:", na.rm, "\n")
  cat("overwrite:", overwrite, "\n")
  
  if(length(scMuffinList$normalized)==0){
    stop("scMuffinList does not contain genes_by_cells\n")
  }
  
  idx_zero <- which(rowSums(scMuffinList$normalized)==0)
  if(length(idx_zero)>0){
    cat("Found genes with all-zero values. These genes may cause issues. Trying to perform the analysis removing these genes from the gene sets.\n")
  }
  
  if(!overwrite & length(scMuffinList$gene_set_scoring)>0){
    shared_columns <- intersect(colnames(scMuffinList$gene_set_scoring$summary), names(gs_list))
    if(length(shared_columns)>0){
      stop("names(gs_list) must be different from colnames(scMuffinList$gene_set_scoring$summary) when overwrite is FALSE\n")
    }
  }
  
  if(kmin>k){
    message("adjusting kmin to k\n")
    kmin <- k
  }
  
  score_type <- score_type[1]
  
  print(names(gs_list))
  cat("# of gene_sets: ", length(gs_list), "\n")
  
  #dataset bins 
  data_bins <- NULL
  if(null_model){
    cat("defining bins...\n")
    data_bins <- sc_data_bin(scMuffinList$normalized, nbins = nbins, na.rm=na.rm)
  }
  
  #signatures-by-cell matrix
  cat("calculating gene set scores...\n")
  if(mc.cores==1){
    
    #debug...
    #res_signatures <- vector("list", length(signatures))
    #for(i in 1:length(signatures)){
    #		cat(i)
    #		res_signatures[[i]] <- gene_set_score(signatures[[i]], genes_by_cells = as.matrix(genes_by_cells@assays$RNA@data), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin)
    #	}
    
    gs_scores <- lapply(gs_list, function(i_marker_set) gs_score(i_marker_set, genes_by_cells = scMuffinList$normalized, bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin, verbose=verbose, na.rm=na.rm))
    
  }else{
    
    gs_scores <- parallel::mclapply(gs_list, function(i_marker_set) gs_score(i_marker_set, genes_by_cells = as.matrix(scMuffinList$normalized), bins = data_bins, k=k, nmark_min = nmark_min, ncells_min = ncells_min, null_model=null_model, kmin = kmin, verbose=verbose, na.rm=na.rm), mc.cores = mc.cores)
    
  }
  
  
  cat("assembling gene_sets_by_cells matrix...")
  #lapply(gene_set_scores, function(x) print(head(x)))
  
  if(score_type == "relative"){
    gss_by_cell_matrix <- t(do.call(cbind, lapply(gs_scores, function(x) array(x$avg_delta_score, dimnames = list(rownames(x))))))
  }else{
    gss_by_cell_matrix <- t(do.call(cbind, lapply(gs_scores, function(x) array(x$case, dimnames = list(rownames(x))))))
  }
  
  if(length(scMuffinList$gene_set_scoring)==0 | overwrite){
    
    scMuffinList$gene_set_scoring <- list(
      summary=as.data.frame(t(gss_by_cell_matrix)),
      full=gs_scores
    )
    
  }else{
    scMuffinList$gene_set_scoring$summary <- merge(scMuffinList$gene_set_scoring$summary, as.data.frame(t(gss_by_cell_matrix)), by=0, all=T, sort=F)
    rownames(scMuffinList$gene_set_scoring$summary) <- scMuffinList$gene_set_scoring$summary$Row.names
    scMuffinList$gene_set_scoring$summary$Row.names <- NULL
    
    scMuffinList$gene_set_scoring$full <- c(scMuffinList$gene_set_scoring$full, gs_scores) ##append 
  }
  
  return(scMuffinList)
  
}
