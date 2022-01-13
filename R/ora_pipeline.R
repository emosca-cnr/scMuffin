#' ora pipeline
#' @param deg_list gene list of interest
#' @param universe all genes
#' @param gs list of gene lists
#' @param gsid2name data.frame for annotation of gene sets. It must contain the column "gsid" with gene set ids that match those in gs
#' @param mc.cores number of cores
#' @param eg2sym data.frame for translating gene ids into symbol. It must contain the columncs gene_id and symbol
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @param out_dir output directory
#' @param write_tables whether to write or not the output
#' @param id project name
#' @importFrom parallel mclapply
#' @export

ora_pipeline <- function(deg_list=NULL, universe=NULL, gs=NULL, gsid2name=NULL, mc.cores=1, eg2sym=NULL, min_size = 5, max_size = 500, out_dir="./", write_tables=FALSE, id="ora"){
  
  #### ORA
  get_genes <- function(gene_set=NULL, wb=NULL, eg2sym=NULL){
    
    ans <- gene_set[gene_set %in% wb]
    if(!is.null(eg2sym)){
      ans <- sort(eg2sym$symbol[eg2sym$gene_id %in% ans])
    }
    
    ans <- paste0(ans, collapse = ";")
    return(ans)
    
  }
  
  #filter gene set lists according to the universe
  gs_i <- lapply(gs, function(x) filter_gsl(x, universe, min_size = min_size, max_size = max_size))
  
  #ORA
  ora_res <- mclapply(gs_i, function(x) ora(deg_list, universe[!universe %in% deg_list], gsl = x, p_adj_method = "fdr"), mc.cores = mc.cores)
  
  #add pathwway info
  ora_res <- lapply(ora_res, function(x) merge(gsid2name, x, by.x="gsid", by.y="id", all.y=T))
  
  #add gene info
  cat("adding gene symbols...\n")
  #k are the pathway dbs
  for(k in 1:length(ora_res)){
    ora_res[[k]]$wbd_symbols <- unlist(parallel::mclapply(ora_res[[k]]$gsid, function(x) get_genes(gs_i[[k]][names(gs_i[[k]]) == x][[1]], deg_list, eg2sym), mc.cores=mc.cores))
  }
  
  #re-order and attach the project id
  ora_res <- lapply(ora_res, function(x) data.frame(id=id, x[order(x$p_adj, x$p, -x$er), ], stringsAsFactors = F))
  
  #write the output
  if(write_tables){
    dir.create(out_dir, recursive = T)
    for(k in 1:length(ora_res)){
      file <- paste0(out_dir, "/", id, "_", names(ora_res)[k], ".txt")
      write.table(ora_res[[k]], file=file, sep="\t", row.names = F)
    }
  }
  
  return(ora_res)
}
