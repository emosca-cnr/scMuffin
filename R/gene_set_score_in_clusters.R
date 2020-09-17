#' gene_set_score_in_clusters
#' @param current_gene_set real data
#' @param seurat_object_data control data
#' @param seurat_object_ident min number of markers
#' @param nmark_min min number of cells
#' @param perc min percentage of gene set
#' @param ncells_min min number of cells
#' @param bins data bins
#' @param k number of permutations
#' @param alt alterative
#' @param test type of test

gene_set_score_in_clusters <- function(current_gene_set, seurat_object_data, seurat_object_ident, nmark_min=5, perc=0.75, ncells_min = NULL, bins = NULL, k=100, alt="g", test="t"){

  clusters <- unique(seurat_object_ident)

  if(is.null(ncells_min)){
    ncells_min <- length(clusters) * 50
  }

  #check whether current gene set is available in input data
  nmark <- sum(rownames(seurat_object_data) %in% current_gene_set)
  curr_perc <- nmark / length(current_gene_set)

  cat("cluster marker availability:", nmark, "/", length(current_gene_set), ";", curr_perc, "\n")
  if(curr_perc > perc & nmark >= nmark_min){

    #null of the current gene set in current subject
    cat("creating null...\n")
    cl1_null <- sc_create_null(seurat_object_data, bins = bins, gene_set = current_gene_set, k = k)

    null_mod_gs <- list(gs=hist(log2(rowSums(seurat_object_data[rownames(seurat_object_data) %in% current_gene_set, ])), 20, plot = F))
    null_mod_gs <- c(null_mod_gs, lapply(cl1_null, function(x) hist(log2(rowSums(x)), 20, plot = F)))

    cat("calculatin gene set score...\n")
    gs_score_res <- gene_set_score(gene_set_data = seurat_object_data[rownames(seurat_object_data) %in% current_gene_set, , drop=F], control_set_data = cl1_null, nmark_min = nmark_min, ncells_min = ncells_min)

    cells_ok <- gs_score_res$cells_ok

    #add cluster information for each cell/barcode
    if(is.data.frame(gs_score_res$res)){

      #delta_scores <- merge(seurat_object_ident, gs_score_res$res, by=0, sort=F)
      delta_scores <- merge(data.frame(seurat_object_ident, stringsAsFactors = F), gs_score_res$res, by=0, sort=F) #modified 2020-06-09; keep all clusters

      colnames(delta_scores)[1:2] <- c("barcode", "cluster")

      #remove clusters with a number of cells less than 10
      cluster_cells <- table(delta_scores$cluster)
      delta_scores <- delta_scores[delta_scores$cluster %in% names(cluster_cells)[cluster_cells >= ncells_min], ]

      if(nrow(delta_scores) >= ncells_min){ #at least 1 cluster with min number of cells

        delta_score_cluster <- data.frame(
          cells=tapply(delta_scores$case, delta_scores$cluster, function(x) length(x)),
          case=tapply(delta_scores$case, delta_scores$cluster, function(x) median(x, na.rm = T)),
          control=tapply(delta_scores$avg_control, delta_scores$cluster, function(x) median(x, na.rm = T)),
          score=tapply(delta_scores$avg_delta_score, delta_scores$cluster, function(x) median(x, na.rm = T)),
          case.z=tapply(delta_scores$case.z, delta_scores$cluster, function(x) median(x, na.rm = T)),
          stringsAsFactors = F
        ) #median cell score for markers(i) in each cell cluster

        if(test == "wrs"){
          temp <- by(delta_scores, delta_scores$cluster, function(x) wilcox.test(x$case, x$avg_control, alternative = alt, paired = T))
        }
        if(test== "t"){
          #temp <- by(delta_scores, delta_scores$cluster, function(x) ifelse(nrow(x) >= 10 & sum(!is.na(x$avg_control)) >= 10, t.test(x$case, x$avg_control, alternative = alt, paired = T), NA))

          temp <- split(delta_scores, delta_scores$cluster)
          for(i in 1:length(temp)){
            if(nrow(temp[[i]]) >= 10 & sum(!is.na(temp[[i]]$avg_control)) >= 10){
              temp[[i]] <- t.test(temp[[i]]$case, temp[[i]]$avg_control, alternative = alt, paired = T)
            }else{
              temp[[i]] <- NA
            }
          }

        }
        res_stat <- unlist(lapply(temp, function(x) ifelse(class(x) == "htest", as.numeric(x$statistic), 0)))
        res_p <- unlist(lapply(temp, function(x) ifelse(class(x) == "htest", as.numeric(x$p.value), 1)))

        delta_score_cluster <- merge(delta_score_cluster, res_stat, by=0)
        delta_score_cluster <- merge(delta_score_cluster, res_p, by.x=1, by.y=0)

        colnames(delta_score_cluster) <- c("cluster", "cells", "case", "control", "score", "case.z", "stat", "p")
        delta_score_cluster$fdr <- p.adjust(delta_score_cluster$p, method = "fdr")

        #missing values to 0
        if(!all(clusters %in% delta_score_cluster$cluster)){
          clusters_missing <- clusters[!clusters %in% delta_score_cluster$cluster]
          delta_score_cluster <- rbind(delta_score_cluster, data.frame(cluster=clusters_missing, cells=NA, case=NA, control=NA, score=0, case.z=NA, stat=0, p=1, fdr=1, stringsAsFactors = F))

        }


      }else{
        delta_score_cluster <- data.frame(cluster=clusters, cells=NA, case=NA, control=NA, score=0, case.z=NA, stat=0, p=1, fdr=1, stringsAsFactors = F)
        delta_scores <- NA

      }


      #for output
      gs_score_res <- gs_score_res$gs_score


    }else{
      delta_score_cluster <- data.frame(cluster=clusters, cells=NA, case=NA, control=NA, score=0, case.z=NA, stat=0, p=1, fdr=1, stringsAsFactors = F)
      gs_score_res <- NA
      delta_scores <- NA
    }

  }else{
    delta_score_cluster <- data.frame(cluster=clusters, cells=NA, case=NA, control=NA, score=0, case.z=NA, stat=0, p=1, fdr=1, stringsAsFactors = F)
    gs_score_res <- NA
    delta_scores <- NA
  }

  #ans <- list(score_clusters=delta_score_cluster, df=delta_scores, perc=curr_perc, nmark=nmark, null_mod_gs=null_mod_gs, gs_score_res=gs_score_res, cells_ok=cells_ok)

  ans <- list(score_clusters=delta_score_cluster, df=delta_scores, perc=curr_perc, nmark=nmark, null_mod_gs=null_mod_gs, gs_score_res=gs_score_res[[1]], cells_ok=cells_ok) #only real

  return(ans)


}
