plot_cluster_Topfeatures<- function(single_cell, cl_lab, output){

qc_data <- FetchData(single_cell, vars = c( "ident", "UMAP_1", "UMAP_2"))
qc_data$cell_ID <- rownames(qc_data)
qc_data_split <- split(qc_data, qc_data$ident)
cl_lab <- cl_lab[match(names(qc_data_split), names(cl_lab))]
topN <- 2
df <- data.frame(cl= length(qc_data_split), cc1 = 0, cc2 = 0, stringName = 0, stringsAsFactors = F)
for(i in 1:length(qc_data_split)){
  temp <- qc_data_split[[i]]
  df[i,] <- c(names(qc_data_split)[i], median(temp$UMAP_1), median(temp$UMAP_2), paste0("\n\n\n",cl_lab[[i]][1], "\n",cl_lab[[i]][2]))
}

png(output, width=200, height=200,  units="mm", res=300)
p <- ggplot(qc_data) + 
  theme_bw()+
  geom_point(aes(x = UMAP_1, y=UMAP_2, color = ident))+
  labs(color="Cluster ID")+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  annotate("text", x = as.numeric(df$cc1), y=as.numeric(df$cc2),  label =df$cl,size=6, fontface=2)+
  annotate("text", x = as.numeric(df$cc1), y=as.numeric(df$cc2),  label =df$stringName,size=2)
#annotate("text", x = as.numeric(df2$cc1), y=as.numeric(df2$cc2)+abs(as.numeric(df2$cc2))/10,  label = df2$V2)
print(p)
dev.off()
}