boxplot_cluster <- function(matrix_features, single_cell, ntop, dir_out){
  
  #create seuratObject, data normalization
  matrix_SO <- CreateSeuratObject(counts = matrix_features, min.cells = 0, min.features = 0)
  all.genes <- rownames(matrix_SO)
  matrix_SO <- ScaleData(matrix_SO, features = all.genes)
  matrix_Scaled <-GetAssayData(object = matrix_SO, slot = "scale.data")
  
  #extract cluster information from single cell data
  qc_data <- FetchData(single_cell, vars = c( "ident", "UMAP_1", "UMAP_2"))
  qc_data$cell_ID <- paste0(rownames(qc_data),"-1")
  qc_data_list <- split(qc_data, qc_data$ident)
  
  #separate cell by cluster ID
  list_cells <- vector("list", length(unique(qc_data$ident)))
  for(i in 1:length(list_cells)){
    list_cells[[i]] <- matrix_Scaled[,which(colnames(matrix_Scaled) %in% qc_data_list[[i]]$cell_ID)]
  }
  names(list_cells) <- names(qc_data_list)
  
  #t-test for each features in each cluster
  mat <- matrix(0, ncol = length(list_cells), nrow = dim(list_cells[[1]])[1])
  for(k in 1:length(list_cells)){
    for(kj in 1:dim(list_cells[[1]])[1]){
      tempList <- lapply(list_cells, function(x) x[kj,])
      feat_Cl <- unlist(tempList[[k]])
      tempList <- tempList[-k]
      feat_notCl <- unlist(tempList)
      mat[kj,k]<-t.test(feat_Cl, feat_notCl)$p.value
    }
  }
  
  
  cl_lab <- matrix(0, nrow = length(list_cells), ncol = ntop+1)
  cl_lab[,1] <- names(list_cells)
  for(cl in 1:length(list_cells)){
    id_lab <- which(mat[,cl] <= sort(mat[,cl], decreasing=F)[ntop], arr.ind=TRUE)
    for(nn in 1:ntop){
      cl_lab[cl,nn+1]  <- rownames(matrix_Scaled)[id_lab[nn]]
    }
  }
  dir_output <- paste0(dir_out,"/boxplot_cluster/")
  dir.create( dir_output )
  for(cl in 1:length(list_cells)){
    png(paste0(dir_output,"clust_",names(list_cells)[cl],".png"), width=380, height=180,  units="mm", res=300)
    par(mar = c(10,4,2,1))
    matrix_Scaled_noCl <- matrix_Scaled[,-which(colnames(matrix_Scaled) %in% colnames(list_cells[[cl]]))]
    boxplot(t(matrix_Scaled_noCl), las=2, main = paste0("Cluster ",names(list_cells)[cl]), xaxt = "n", pch=16)
    temp_cells <- list_cells[[cl]]
    index <- which(mat[,cl] <= sort(mat[,cl], decreasing=F)[5], arr.ind=TRUE)
    index_temp <- which(mat[,cl] <= 0.01, arr.ind=TRUE)
    index2 <- setdiff(index_temp,index)
    index3 <- setdiff(seq(1,dim(mat)[1],1), index_temp )
    setdiff(index2,index)
    for(ji in 1:length(index)){
      points(jitter(rep(index[ji],dim(temp_cells)[2]), amount=0.3), temp_cells[index[ji],], col=adjustcolor("red",0.7), pch=16, cex=0.6)
    }
    for(ji in 1:length(index2)){
      points(jitter(rep(index2[ji],dim(temp_cells)[2]), amount=0.3), temp_cells[index2[ji],], col=adjustcolor("salmon",0.3), pch=16, cex=0.6)
    } 
    for(ji in 1:length(index3)){
      points(jitter(rep(index3[ji],dim(temp_cells)[2]), amount=0.3), temp_cells[index3[ji],], col=adjustcolor("peachpuff",0.1), pch=16, cex=0.6)
    }   
    text(seq(1:dim(temp_cells)[1]), par("usr")[3] - 0.3, labels = rownames(matrix_Scaled), srt = 45, adj = 1, xpd = TRUE)
    dev.off()
  }
}
