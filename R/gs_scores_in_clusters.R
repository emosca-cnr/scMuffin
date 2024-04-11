#' Gene set scoring at cluster level
#' @param score_table real data
#' @param cell_clusters vector of cell clusters
#' @param ncells_min minimum number of cells required for the calculation of the average signature in the cluster
#' @param alt alterative passed to [wilcox.test()] or [t.test()]
#' @param test type of test: t to use [t.test()]; wrs to use [wilcox.test()]
#' @param fract_min only clusters with this fraction of cells with not null gene set score will be considered
#' @description Gene set scoring in clusters
#' @importFrom stats median wilcox.test t.test p.adjust
#' @export

gs_scores_in_clusters <- function(score_table=NULL, cell_clusters=NULL, ncells_min=5, alt="g", test="t", fract_min=0.5){
	
	if(!is.factor(cell_clusters)){
		cell_clusters <- as.factor(cell_clusters)
	}
	
	null_model <- FALSE
	if(any(!is.na(score_table$avg_control))){
		null_model <- TRUE
	}
	
  clusters_size <- table(cell_clusters)
	clusters <- levels(cell_clusters)
	
	score_table_clusters <- merge(data.frame(cluster=as.character(cell_clusters), stringsAsFactors = F, row.names = names(cell_clusters)), score_table, by=0, sort=F) #modified 2020-06-09; keep all clusters
	colnames(score_table_clusters)[1:2] <- c("cell", "cluster")
	
	if(null_model){
		score_table_clusters$filter <- score_table_clusters$nmark_min & score_table_clusters$null_ok
	}else{
		score_table_clusters$filter <- score_table_clusters$nmark_min
	}
	score_table_clusters <- score_table_clusters[score_table_clusters$filter, ]
	
	if(nrow(score_table_clusters)>= ncells_min){ 
		
		if(null_model){
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
			
			cluster_scores <- data.frame(
				cells=tapply(score_table_clusters$nmark_min, score_table_clusters$cluster, sum),
				med_case=tapply(score_table_clusters$case, score_table_clusters$cluster, function(x) median(x, na.rm = T)),
				med_control=NA,
				stringsAsFactors = F
			) #median cell score for markers(i) in each cell cluster
			
			cluster_scores$score <- cluster_scores$med_case
			
			cluster_scores$stat <- 0
			cluster_scores$p1 <- 1
			cluster_scores$fdr <- 1
			
			cluster_scores <- data.frame(cluster=rownames(cluster_scores), cluster_scores, stringsAsFactors = F)
			
		}
		
	}else{
		warning("not enough cells in any cluster")
		cluster_scores <- data.frame(cluster=clusters, cells=NA, med_case=NA, med_control=NA, score=0, stat=0, p=1, fdr=1, stringsAsFactors = F)
	}
	
	#cluster original size
	cluster_scores$size <- clusters_size[match(cluster_scores$cluster, names(clusters_size))]
	
	#requirements
	idx_na <- ((cluster_scores$cells / cluster_scores$size) < fract_min) | (cluster_scores$cells < ncells_min)
	cluster_scores$med_case[idx_na] <- NA
	cluster_scores$med_control[idx_na] <- NA
	cluster_scores$score[idx_na] <- 0
	cluster_scores$stat[idx_na] <- 0
	cluster_scores$p[idx_na] <- 0
	cluster_scores$fdr[idx_na] <- 0
	
	cluster_scores <- cluster_scores[match(clusters, cluster_scores$cluster), ]
	
	return(cluster_scores)
	
}
