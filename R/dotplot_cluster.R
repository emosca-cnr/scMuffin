#' Boxplot clusters 
#' 
#' Produce boxplots of the given features in each cluster. A t-test is performed for each feature among clusters.
#' @param features feautures object
#' @param clusterings clusterings object
#' @param clustering_name one among the clusterings
#' @param clust_enrich_res clustering enrichment results
#' @param text_val which column of the hypergeometric result. Default is "p". Other values are "wbd" and "er"
#' @param dir_out string, output directory
#' @importFrom grDevices jpeg
#' @export

dotplot_cluster <- function(features, clusterings=NULL, clustering_name=NULL, dir_out="./", clust_enrich_res=NULL, text_val=c("p", "wbd", "er")){
	
  if(!dir.exists(dir_out)){
    dir.create(dir_out, recursive = TRUE)
  }
  
  text_val <- match.arg(text_val, c("p", "wbd"))
  
  cells_by_features <- features$df[, features$type=="factor", drop=F]
  
  cell_clusters <- clusterings[, colnames(clusterings) == clustering_name]
  cell_clusters_set <- levels(cell_clusters)
	
  # clusters-by-feature value table of enrichment p-value
  en_res <- extract_cluster_enrichment_table(clust_enrich_res, c_type = text_val)
  en_res <- en_res$cluster_hyper_table[names(en_res$cluster_hyper_table) == clustering_name][[1]]
  
  if(text_val == "wbd"){
    for(i in 1:length(en_res)){
      en_res[[i]] <- t(apply(en_res[[i]], 1, function(x) x/sum(x)))
    }
  }
  
	#boxplot for each cluster
	for(cl in 1:length(cell_clusters_set)){ ###for each cluster
		
		grDevices::jpeg(paste0(dir_out, "/cluster_", cell_clusters_set[cl],".jpg"), width=180, height=180, units="mm", res=300)
		par(mar = c(10, 4, 2, 1))
		
		#distribution of all cells by feature
		#feature data of the cluster
		data_clust <- as.data.frame(cells_by_features[rownames(cells_by_features) %in% rownames(clusterings)[cell_clusters == cell_clusters_set[cl]], ])
		
		#feature data other clusters
		data_no_clust <- as.data.frame(cells_by_features[!rownames(cells_by_features) %in% rownames(clusterings)[cell_clusters == cell_clusters_set[cl]], ])
		
		#boxplots of the cluster
		data_clust_at <- seq(1, ncol(cells_by_features)*2, by = 2)
		data_no_clust_at <- seq(2, ncol(cells_by_features)*2, by = 2)
		
		plot(0, xlim=c(0, max(data_no_clust_at))+0.5, pch="", ylim = c(1, max(sapply(cells_by_features, as.numeric))), xaxt="n", xlab="", ylab="", main = paste0("Cluster ", cell_clusters_set[cl]), yaxt="n")
		
		for(i in 1:ncol(cells_by_features)){ #for each feature
			
		#   col <- "pink"
		# 	if(!is.null(top_features)){
		# 		if(colnames(cells_by_features)[i] %in% gsub("_.+$", "", unlist(top_features[names(top_features) == cell_clusters_set[cl]]))){
		# 			col <- "red"
		# 		}
		# 	}
			col <- "red"
			
			data_clust_i_fact <- data_clust[, i]
			data_clust_i <- as.numeric(data_clust_i_fact)
			
			points(jitter(rep(data_clust_at[i], nrow(data_clust)), amount=0.3), jitter(data_clust_i, amount = 0.3), col=adjustcolor(col, 0.7), pch=16, cex=0.3)
			j <- (i-1)*length(cell_clusters_set)+cl
			cat(j)
			
			#obs_f <- cont_tables[[j]]$observed / sum(cont_tables[[j]]$observed)
			#obs_f <- obs_f[order(as.numeric(names(obs_f)))]
			#text(rep(data_clust_at[i], length(obs_f)), as.numeric(as.factor(names(obs_f))), format(obs_f, digits=2), cex=0.7, font=2)
			#text(rep(data_clust_at[i]+0.5, length(obs_f)), as.numeric(as.factor(names(obs_f))), names(obs_f), cex=0.7)
			
			p_val_cl_i <- en_res[[i]][rownames(en_res[[i]]) == cell_clusters_set[cl], ]
			text(rep(data_clust_at[i], length(p_val_cl_i)), as.numeric(as.factor(names(p_val_cl_i))), format(p_val_cl_i, digits=2), cex=0.7, font=2)
			
			text(rep(data_clust_at[i]+0.5, length(p_val_cl_i)), 1:length(levels(data_clust_i_fact)), levels(data_clust_i_fact), cex=0.7, font=2, col="red")
			
			
		}
		
		
		#dotplots of other clusters
		for(i in 1:ncol(cells_by_features)){
			
			points(jitter(rep(data_no_clust_at[i], nrow(data_no_clust)), amount=0.3), jitter(as.numeric(data_no_clust[, i])), col=adjustcolor("darkgray", 0.7), pch=16, cex=0.3)
			j <- (i-1)*length(cell_clusters_set)+cl
			cat(j)
			
		}
		
		axis(1, data_clust_at+0.5, colnames(cells_by_features), las=2, cex=0.5)
		
		dev.off()
	}
}
