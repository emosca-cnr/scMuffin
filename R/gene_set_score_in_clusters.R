#' gene_set_score_in_clusters
#' @param current_gene_set real data
#' @param seurat_object_data control data
#' @param seurat_object_ident min number of markers
#' @param alt alterative
#' @param test type of test
#' @import stats
#' @export

gene_set_score_in_clusters <- function(score_table, cell_clusters, ncells_min=5, alt="g", test="t"){
	
	clusters <- unique(cell_clusters)
	
	score_table_clusters <- merge(data.frame(cluster=cell_clusters, stringsAsFactors = F), score_table, by=0, sort=F) #modified 2020-06-09; keep all clusters
	colnames(score_table_clusters)[1:2] <- c("cell", "cluster")
	
	score_table_clusters$filter <- score_table_clusters$nmark_min & score_table_clusters$null_ok
	score_table_clusters <- score_table_clusters[score_table_clusters$filter, ]
	
	if(nrow(score_table_clusters)>= ncells_min){ 
		

	cluster_scores <- data.frame(
		cells=tapply(score_table_clusters$nmark_min, score_table_clusters$cluster, sum),
		med_case=tapply(score_table_clusters$case, score_table_clusters$cluster, function(x) median(x, na.rm = T)),
		med_control=tapply(score_table_clusters$avg_control, score_table_clusters$cluster, function(x) median(x, na.rm = T)),
		stringsAsFactors = F
	) #median cell score for markers(i) in each cell cluster
	
	cluster_scores$score <- cluster_scores$med_case - cluster_scores$med_control
	
	if(test == "wrs"){
		
		temp <- split(score_table_clusters, score_table_clusters$cluster)
		for(i in 1:length(temp)){
			if(nrow(temp[[i]]) >= ncells_min ){ #at least ncells_min with acceptable score
				temp[[i]] <- wilcox.test(temp[[i]]$case, temp[[i]]$avg_control, alternative = alt, paired = F)
			}else{
				temp[[i]] <- NA
			}
		}
	}
	
	if(test== "t"){
		
		temp <- split(score_table_clusters, score_table_clusters$cluster)
		for(i in 1:length(temp)){
			if(nrow(temp[[i]]) >= ncells_min ){ #at least ncells_min with acceptable score
				temp[[i]] <- t.test(temp[[i]]$case, temp[[i]]$avg_control, alternative = alt, paired = T)
			}else{
				temp[[i]] <- NA
			}
		}
		
	}
	res_stat <- unlist(lapply(temp, function(x) ifelse(class(x) == "htest", as.numeric(x$statistic), 0)))
	res_p <- unlist(lapply(temp, function(x) ifelse(class(x) == "htest", as.numeric(x$p.value), 1)))
	
	cluster_scores <- merge(cluster_scores, res_stat, by=0)
	cluster_scores <- merge(cluster_scores, res_p, by.x=1, by.y=0)
	
	colnames(cluster_scores)[c(1, 6, 7)] <- c("cluster", "stat", "p")
	cluster_scores$fdr <- p.adjust(cluster_scores$p, method = "fdr")
	
	#missing values to 0
	if(!all(clusters %in% cluster_scores$cluster)){
		clusters_missing <- clusters[!clusters %in% cluster_scores$cluster]
		cluster_scores <- rbind(cluster_scores, data.frame(cluster=clusters_missing, cells=NA, med_case=NA, med_control=NA, score=0, stat=0, p=1, fdr=1, stringsAsFactors = F))
		
	}
	
	}else{
		warning("not enough cells in any cluster")
		cluster_scores <- data.frame(cluster=clusters, cells=NA, med_case=NA, med_control=NA, score=0, stat=0, p=1, fdr=1, stringsAsFactors = F)
	}
	
	return(cluster_scores)

}
