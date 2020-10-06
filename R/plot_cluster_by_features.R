#' plot_cluster_by_features
#'
#'
#'

plot_cluster_by_features <- function(features_by_cells, matrix_Scaled, features_by_cells, cl_lab, output){
	
	t_matrix_Scaled <- t(matrix_Scaled)
	
	#colnames(t_matrix_Scaled) <- paste0("`", colnames(t_matrix_Scaled), "`")
	qc_data <- FetchData(features_by_cells, vars = c( "ident", "UMAP_1", "UMAP_2"))
	qc_data$cell_ID <- paste0(rownames(qc_data),"-1")
	cell_features <- merge(qc_data, t_matrix_Scaled,by.x="cell_ID", by.y=0)
	
	qc_data_split <- split(qc_data, qc_data$ident)
	cl_lab <- cl_lab[match(names(qc_data_split), names(cl_lab))]
	df <- data.frame(cl= length(qc_data_split), cc1 = 0, cc2 = 0, stringName = 0, stringsAsFactors = F)
	for(i in 1:length(qc_data_split)){
		temp <- qc_data_split[[i]]
		df[i,] <- c(names(qc_data_split)[i], median(temp$UMAP_1), median(temp$UMAP_2), paste0("\n\n\n",cl_lab[[i]][1], "\n",cl_lab[[i]][2]))
	}
	
	dir_output <- paste0(output,"/cluster_byFeatures/")
	dir.create( dir_output )
	for(j in 5:dim(cell_features)[2]){
		png(paste0(dir_output,"cluster_by_",colnames(cell_features)[j],".png"),  width=200, height=200,  units="mm", res=300)
		var_n <- "`colnames(cell_features)[j]`"
		p<- ggplot(cell_features, aes_string(x="UMAP_1", y="UMAP_2", color=paste0("`",colnames(cell_features)[j],"`"))) + 
			geom_point()+ggtitle(paste0("Colored by ",colnames(cell_features)[j]))+
			annotate("text", x = as.numeric(df$cc1), y=as.numeric(df$cc2),  label =df$cl,size=6, fontface=2, col="black")+
			annotate("text", x = as.numeric(df$cc1), y=as.numeric(df$cc2),  label =df$stringName,size=2, col="black")
		pcol <- scale_color_gradientn(colors =  c("lightgrey", "blue"))
		q <-  p+pcol +theme_bw()
		print(q)
		dev.off()
	}
}
