#' #' # # find_variable_features <- function (object, selection.method = "vst", loess.span = 0.3, clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR, num.bin = 20, binning.method = "equal_width", verbose = TRUE){
#' #' # # 	
#' #' # # 	clip.max <- sqrt(x = ncol(x = object)) #sqrt of the number of cells
#' #' # # 	
#' #' # # 	hvf.info <- data.frame(mean = rowMeans(x = object))
#' #' # # 	hvf.info$variance <- SparseRowVar2(mat = object, mu = hvf.info$mean, display_progress = verbose)
#' #' # # 	hvf.info$variance.expected <- 0
#' #' # # 	hvf.info$variance.standardized <- 0
#' #' # # 	not.const <- hvf.info$variance > 0
#' #' # # 	fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
#' #' # # 							 data = hvf.info[not.const, ], span = loess.span)
#' #' # # 	hvf.info$variance.expected[not.const] <- 10^fit$fitted
#' #' # # 	hvf.info$variance.standardized <- SparseRowVarStd(mat = object, 
#' #' # # 																										mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected), 
#' #' # # 																										vmax = clip.max, display_progress = verbose)
#' #' # # 	colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
#' #' # # }
#' #' # # else {
#' #' # # 	if (!inherits(x = mean.function, what = "function")) {
#' #' # # 		stop("'mean.function' must be a function")
#' #' # # 	}
#' #' # # 	if (!inherits(x = dispersion.function, what = "function")) {
#' #' # # 		stop("'dispersion.function' must be a function")
#' #' # # 	}
#' #' # # 	feature.mean <- mean.function(object, verbose)
#' #' # # 	feature.dispersion <- dispersion.function(object, verbose)
#' #' # # 	names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
#' #' # # 	feature.dispersion[is.na(x = feature.dispersion)] <- 0
#' #' # # 	feature.mean[is.na(x = feature.mean)] <- 0
#' #' # # 	data.x.breaks <- switch(EXPR = binning.method, equal_width = num.bin, 
#' #' # # 													equal_frequency = c(-1, quantile(x = feature.mean[feature.mean > 
#' #' # # 																																							0], probs = seq.int(from = 0, to = 1, length.out = num.bin))), 
#' #' # # 													stop("Unknown binning method: ", binning.method))
#' #' # # 	data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
#' #' # # 	names(x = data.x.bin) <- names(x = feature.mean)
#' #' # # 	mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, 
#' #' # # 									 FUN = mean)
#' #' # # 	sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, 
#' #' # # 								 FUN = sd)
#' #' # # 	feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)])/sd.y[as.numeric(x = data.x.bin)]
#' #' # # 	names(x = feature.dispersion.scaled) <- names(x = feature.mean)
#' #' # # 	hvf.info <- data.frame(feature.mean, feature.dispersion, 
#' #' # # 												 feature.dispersion.scaled)
#' #' # # 	rownames(x = hvf.info) <- rownames(x = object)
#' #' # # 	colnames(x = hvf.info) <- paste0("mvp.", c("mean", "dispersion", 
#' #' # # 																						 "dispersion.scaled"))
#' #' # # }
#' #' # # return(hvf.info)
#' #' # # }
#' #' # # 
#' #' # # if (selection.method == "vst") {
#' #' # # 	hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, 
#' #' # # 														 decreasing = TRUE), , drop = FALSE]
#' #' # # }
#' #' # 
#' #' # 
#' #' # find_variable_features <- function (features_by_cells, nfeatures = 50){
#' #' # 	
#' #' # 	ans <- features_by_cells
#' #' # 	if(nrow(ans) <= nfeatures){
#' #' # 		
#' #' # 		warnings("input data lower or equal to ", nfeatures, "\n")
#' #' # 		
#' #' # 	}else{
#' #' # 		
#' #' # 		#scale data
#' #' # 		ans <- CreateSeuratObject(ans, min.cells = 0, min.features = 0)
#' #' # 		ans <- FindVariableFeatures(ans, selection.method = "vst", nfeatures = nfeatures)
#' #' # 		
#' #' # 		
#' #' # 	}
#' #' # 	
#' #' # 	return(ans)
#' #' # 	
#' #' # }
#' #' # 
#' #' 
#' #' 
#' #' # #' heatmap of signatures
#' #' # #' @param X ...
#' #' # #' @param file ...
#' #' # #' @description Plot as an heatmap the resulting data from CNV. 
#' #' # #' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' #' # #' @usage heatmap_CNV(chr_merged)
#' #' # #' @author Ettore Mosca
#' #' # #' @import RColorBrewer graphics
#' #' # ## @importFrom dendextend color_branches
#' #' # #' @export
#' #' # 
#' #' # heatmap_cluster_enrichment_ <- function(X, Y, file="heatmap_cluster_enrichment.jpg", pal=NULL, n_colors=11, seurat_dendrogram=NULL, width=180, height=180, res=300, cex.axis=0.3, fdr_cutoff=0.05) {
#' #' # 	
#' #' # 	rotate <- function(x) t(apply(x, 2, rev)) # rotate +90
#' #' # 	
#' #' # 	if(is.null(pal)){
#' #' # 		colors_ <- rev(brewer.pal(n_colors, "RdYlBu"))
#' #' # 	}
#' #' # 	
#' #' # 	
#' #' # 	#cells clustering
#' #' # 	if(!is.null(seurat_dendrogram)){
#' #' # 		hc_col <- seurat_dendrogram
#' #' # 		X <- X[, match(hc_col$labels, gsub("^.+_", "", colnames(X)))]
#' #' # 	}else{
#' #' # 		hc_col <- hclust(dist(t(X))) #COLUMNS
#' #' # 		X <- X[, hc_col$order]
#' #' # 	}
#' #' # 	
#' #' # 	#signatures clustering
#' #' # 	hc_row <- hclust(dist(X)) #ROWS
#' #' # 	X <- X[hc_row$order, ]
#' #' # 	
#' #' # 	
#' #' # 	sample_labels <- apply(X, 2, which.max)
#' #' # 	
#' #' # 	#clust_col <- cutree(hc_col, 7)
#' #' # 	#clust_col <- clust_col[match(colnames(X), names(clust_col))]
#' #' # 	#clust_col_sep <- which(clust_col[2:length(clust_col)] - clust_col[1:(length(clust_col)-1)]!=0)
#' #' # 	#dend_col <- dendextend::color_branches(as.dendrogram(hc_col), k = 7)
#' #' # 	dend_col <- as.dendrogram(hc_col)
#' #' # 	
#' #' # 	#clust_row <- cutree(hc_row, 5)
#' #' # 	#clust_row <- clust_row[match(rownames(X), names(clust_row))]
#' #' # 	#clust_row_sep <- which(clust_row[2:length(clust_row)] - clust_row[1:(length(clust_row)-1)]!=0)
#' #' # 	#dend_row <- dendextend::color_branches(as.dendrogram(hc_row), k = 5)
#' #' # 	dend_row <- as.dendrogram(hc_row)
#' #' # 	
#' #' # 	#rotation
#' #' # 	X <- rotate(X)
#' #' # 	Y <- rotate(Y)
#' #' # 	Y <- Y[match(rownames(X), rownames(Y)), match(colnames(X), colnames(Y))]
#' #' # 	
#' #' # 	jpeg(file, width = width, height = height, res=res, units="mm")
#' #' # 	
#' #' # 	layout.show(layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = T), heights = c(0.15, 0.95), widths = c(0.2, 0.8)))
#' #' # 	
#' #' # 	par(mar=c(0, 0, 0, 0))
#' #' # 	plot.new()
#' #' # 	par(xpd=T)
#' #' # 	
#' #' # 	leg_text <- cbind(c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-12], c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-1])
#' #' # 	leg_text <- apply(leg_text, 1, function(x) paste0(format(x[1], digits = 2), ", ", format(x[2], digits = 2)))
#' #' # 	
#' #' # 	legend("bottom", leg_text, pch=16, col=colors_, cex=0.4, xpd=T, pt.cex=0.6)
#' #' # 	
#' #' # 	c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))
#' #' # 	
#' #' # 	par(mar=c(0, 2, 1, 0))
#' #' # 	plot(dend_col, type = "rectangle", leaflab = "none", axes=F, edgePar = list())
#' #' # 	
#' #' # 	par(mar=c(1.5, 0, 0, 2))
#' #' # 	plot(rev(dend_row), type = "rectangle", leaflab = "none", axes=F, horiz = T)
#' #' # 	#plot.new()
#' #' # 	
#' #' # 	par(mar=c(6, 3, 1.5, 1))
#' #' # 	image(X, xaxt="none", yaxt="none", col=colors_, breaks = c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6)))
#' #' # 	
#' #' # 	xx <- seq(0, 1, length.out = nrow(X))
#' #' # 	yy <- seq(0, 1, length.out = ncol(X))
#' #' # 	xx_lab <- rownames(X) #horizontal axis
#' #' # 	yy_lab <- colnames(X) #vertical axis
#' #' # 	axis(1, xx, xx_lab, las=2, cex.axis=cex.axis)
#' #' # 	axis(2, yy, yy_lab, las=2, cex.axis=cex.axis)
#' #' # 	
#' #' # 	for(i in 1:nrow(X)){ #xx
#' #' # 		for(j in 1:ncol(X)){ #yy
#' #' # 			
#' #' # 			if(Y[i, j] < fdr_cutoff & X[i, j] > 0){ #significant and positive nes
#' #' # 				points(xx[i], yy[j], pch="*", cex=2, lwd=2, col="pink")
#' #' # 			}
#' #' # 			
#' #' # 		}
#' #' # 	}
#' #' # 	
#' #' # 	dev.off()
#' #' # 	
#' #' # }
#' #' #
#' #' 
#' #' ##' heatmap of signatures
#' #' ##' @param X ...
#' #' ##' @param file ...
#' #' ##' @description Plot as an heatmap the resulting data from CNV. 
#' #' ##' @details Preprocessing with 'preprocess_for_heatmap' needed. 
#' #' ##' @usage heatmap_CNV(chr_merged)
#' #' ##' @author Ettore Mosca
#' #' ##' @import RColorBrewer graphics
#' #' ## @importFrom dendextend color_branches
#' #' ##' @export
#' #' 
#' #' # heatmap_signatures_ <- function(X, file="heatmap_signatures.jpg", pal=NULL, n_colors=11, seurat_dendrogram=NULL, width=180, height=180, res=300, cex.axis=0.3, scale_rows=FALSE) {
#' #' # 	
#' #' # 	rotate <- function(x) t(apply(x, 2, rev)) # rotate +90
#' #' # 	
#' #' # 	if(is.null(pal)){
#' #' # 		colors_ <- rev(brewer.pal(n_colors, "RdYlBu"))
#' #' # 	}
#' #' # 	
#' #' # 	
#' #' # 	if(scale_rows){
#' #' # 		temp <- t(apply(X, 1, scale))
#' #' # 		rownames(temp) <- rownames(X)
#' #' # 		colnames(temp) <- colnames(X)
#' #' # 		X <- temp
#' #' # 		rm(temp)
#' #' # 		X[is.nan(X)] <- 0
#' #' # 	}
#' #' # 	
#' #' # 	#cells clustering
#' #' # 	if(!is.null(seurat_dendrogram)){
#' #' # 		hc_col <- seurat_dendrogram
#' #' # 		X <- X[, match(hc_col$labels, gsub("^.+_", "", colnames(X)))]
#' #' # 	}else{
#' #' # 		hc_col <- hclust(dist(t(X))) #COLUMNS
#' #' # 		X <- X[, hc_col$order]
#' #' # 	}
#' #' # 	
#' #' # 	#signatures clustering
#' #' # 	hc_row <- hclust(dist(X)) #ROWS
#' #' # 	X <- X[hc_row$order, ]
#' #' # 	
#' #' # 	
#' #' # 	sample_labels <- apply(X, 2, which.max)
#' #' # 	
#' #' # 	#clust_col <- cutree(hc_col, 7)
#' #' # 	#clust_col <- clust_col[match(colnames(X), names(clust_col))]
#' #' # 	#clust_col_sep <- which(clust_col[2:length(clust_col)] - clust_col[1:(length(clust_col)-1)]!=0)
#' #' # 	#dend_col <- dendextend::color_branches(as.dendrogram(hc_col), k = 7)
#' #' # 	dend_col <- as.dendrogram(hc_col)
#' #' # 	
#' #' # 	#clust_row <- cutree(hc_row, 5)
#' #' # 	#clust_row <- clust_row[match(rownames(X), names(clust_row))]
#' #' # 	#clust_row_sep <- which(clust_row[2:length(clust_row)] - clust_row[1:(length(clust_row)-1)]!=0)
#' #' # 	#dend_row <- dendextend::color_branches(as.dendrogram(hc_row), k = 5)
#' #' # 	dend_row <- as.dendrogram(hc_row)
#' #' # 	
#' #' # 	#rotation
#' #' # 	X <- rotate(X)
#' #' # 	
#' #' # 	jpeg(file, width = width, height = height, res=res, units="mm")
#' #' # 	
#' #' # 	layout.show(layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = T), heights = c(0.15, 0.95), widths = c(0.2, 0.8)))
#' #' # 	
#' #' # 	par(mar=c(0, 0, 0, 0))
#' #' # 	plot.new()
#' #' # 	par(xpd=T)
#' #' # 	
#' #' # 	leg_text <- cbind(c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-12], c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))[-1])
#' #' # 	leg_text <- apply(leg_text, 1, function(x) paste0(format(x[1], digits = 2), ", ", format(x[2], digits = 2)))
#' #' # 	
#' #' # 	legend("bottom", leg_text, pch=16, col=colors_, cex=0.4, xpd=T, pt.cex=0.6)
#' #' # 	
#' #' # 	c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6))
#' #' # 	
#' #' # 	par(mar=c(0, 2, 1, 0))
#' #' # 	plot(dend_col, type = "rectangle", leaflab = "none", axes=F, edgePar = list())
#' #' # 	
#' #' # 	par(mar=c(1.5, 0, 0, 2))
#' #' # 	plot(rev(dend_row), type = "rectangle", leaflab = "none", axes=F, horiz = T)
#' #' # 	#plot.new()
#' #' # 	
#' #' # 	par(mar=c(3, 3, 1.5, 1))
#' #' # 	image(X, xaxt="none", yaxt="none", col=colors_, breaks = c(seq(min(X), -0.1, length.out = 5), 0, seq(0.1, max(X), length.out = 6)))
#' #' # 	
#' #' # 	xx <- seq(0, 1, length.out = nrow(X))
#' #' # 	yy <- seq(0, 1, length.out = ncol(X))
#' #' # 	xx_lab <- rownames(X) #horizontal axis
#' #' # 	yy_lab <- colnames(X) #vertical axis
#' #' # 	axis(1, xx, xx_lab, las=2, cex.axis=cex.axis)
#' #' # 	axis(2, yy, yy_lab, las=2, cex.axis=cex.axis)
#' #' # 	
#' #' # 	for(i in 1:nrow(X)){ #xx
#' #' # 		for(j in 1:ncol(X)){ #yy
#' #' # 			
#' #' # 			if(X[i, j] == 0){
#' #' # 				text(xx[i], yy[j], "X", cex=.6)
#' #' # 			}
#' #' # 			
#' #' # 			if((ncol(X):1)[sample_labels[i]]==j){
#' #' # 				points(xx[i], yy[j], pch="*", cex=2, lwd=2, col="pink")
#' #' # 			}
#' #' # 		}
#' #' # 	}
#' #' # 	
#' #' # 	dev.off()
#' #' # 	
#' #' # }
#' #' 
#' #' #' preprocess_for_heatmap
#' #' #' @param result_cnv ...
#' #' #' @description Bind chromosomes together as a matrix. 
#' #' #' @usage preprocess_for_heatmap(result_cnv)
#' #' #' @return A matrix to be used by the function heatmap_CNV.
#' #' #' @author Valentina Nale
#' #' #' @export
#' #' 
#' #' preprocess_for_heatmap_CNV <- function(result_cnv) {
#' #' 	
#' #' 	for(i in 1:length(result_cnv)){
#' #' 		rownames(result_cnv[[i]]) <- paste0("chr", names(result_cnv)[i], "_", rownames(result_cnv[[i]]))
#' #' 	} 
#' #' 	
#' #' 	chr_merged <- do.call(rbind, result_cnv)
#' #' 	
#' #' 	return(chr_merged)
#' #' }
#' #' 
#' 
#' 
#' #' preprocess_object_for_CNV
#' #' @param input_object genes-by-cells input matrix
#' #' @description Preprocessing function to obtain a genomically-ordered list of chromosomes. 
#' #' @details The preliminary step consist of annotation, duplicates and NA values removal. 
#' #' Then, the matrix is splitted as a list of dataframe, where every dataframe is a chromosome.
#' #' Chromosomes are ordered from 1 to 22 + X +Y, and then re-ordered by start position. 
#' #' @usage preprocess_object_for_cnv(input_object)
#' #' @return list of genomically-ordered chromosomes
#' #' @author Valentina Nale
#' #' @import Seurat org.Hs.eg.db
#' #' @export
#' preprocess_object_for_CNV <- function(input_object) {
#' 	
#' 	# retrieve gene informations
#' 	# input_object=cellObj
#' 	
#' 	cat("Retrieving gene locations...\n")
#' 	gene_locations <- as.data.frame(org.Hs.eg.db::org.Hs.egCHRLOC)
#' 	temp <- as.data.frame(org.Hs.eg.db::org.Hs.egCHRLOCEND)
#' 	eg2sym <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
#' 	
#' 	gene_locations <- merge(eg2sym, gene_locations, by="gene_id", sort=F)
#' 	gene_locations <- merge(gene_locations, temp, by=c("gene_id", "Chromosome"), sort=F)
#' 	gene_locations$pos <- apply(abs(gene_locations[, c("start_location", "end_location")]), 1, min) 
#' 	gene_locations$start_location <- NULL
#' 	gene_locations$end_location <- NULL
#' 	
#' 	# merged Seurat matrix and annotation matrix via "symbol"
#' 	#input_object <- Seurat::GetAssayData(object=input_object, slot="data")
#' 	matrix_complete <- merge(gene_locations, input_object, by.x="symbol", by.y=0, sort=FALSE)
#' 	
#' 	# remove chromosome names not in 1:22 and X,Y -- long names, duplicates
#' 	matrix_complete <- matrix_complete[-which(grepl("_", matrix_complete$Chromosome)), ]
#' 	
#' 	# remove NA values across entrezid
#' 	matrix_complete <- matrix_complete[!is.na(matrix_complete$gene_id), ]
#' 	
#' 	# remove duplicates (entregene_id and on row.names('symbol'))
#' 	matrix_reduced <- cbind(matrix_complete, mean = rowMeans(matrix_complete[, 5:dim(matrix_complete)[2]]))
#' 	matrix_reduced <- matrix_reduced[order(-matrix_reduced$mean), ]
#' 	matrix_reduced <- matrix_reduced[!duplicated(matrix_reduced$gene_id),]
#' 	matrix_reduced <- matrix_reduced[!duplicated(matrix_reduced$symbol),]
#' 	matrix_reduced$mean <- NULL
#' 	
#' 	# change rownames with 'SYMBOL'
#' 	mat_rn <- matrix_reduced[,-1]
#' 	rownames(mat_rn) <- matrix_reduced[,1]
#' 	
#' 	# split chromosomes into a list of dataframe
#' 	mat_splitted <- split(mat_rn, mat_rn$Chromosome)
#' 	
#' 	# sort chromosomes by position
#' 	mat_sorted <- lapply(mat_splitted, function(df){
#' 		df[order(df$pos),]
#' 	})
#' 	return(mat_sorted)
#' }
#' 
#' 
#' #' #' plot_umap
#' #' 
#' #' Generate a UMAP visualization
#' #' @param seurat_object seurat object, object with saved dimension reduction components
#' #' @param file string, file name output
#' #' @param color_by string, specification of a feature to colour by (e.g. cluster ID)
#' #' 
#' #' @import Seurat graphics
#' #' @export
#' #' 
#' plot_umap_ <- function(seurat_object, file="umap.jpg", color_by="ident", pal=NULL, labels=NULL, ...){
#' 	
#' 	data_plot <- Seurat::FetchData(seurat_object, vars = c("UMAP_1", "UMAP_2", color_by))
#' 	
#' 	col_levels <- levels(factor(data_plot[, 3], levels=sort(unique(data_plot[, 3]))))
#' 	
#' 	jpeg(file, width=180, height=180, units="mm", res=300)
#' 	
#' 	par(mar=c(3, 3, 3, 1))
#' 	par(mgp=c(2, 0.7, 0))
#' 	
#' 	layout(matrix(c(1, 2), nrow = 1), widths = c(0.85, 0.15))
#' 	#plot(data_plot$UMAP_1, data_plot$UMAP_2, pch=16, col=rainbow(length(col_levels))[as.numeric(data_plot[, 3])], cex=0.4, xlab="UMAP1", ylab = "UMAP2")
#' 	plot(data_plot$UMAP_1, data_plot$UMAP_2, pch=16, col=pal[as.numeric(data_plot[, 3])], xlab="UMAP1", ylab = "UMAP2", ...)
#' 	
#' 	if(!is.null(labels)){
#' 		labels <- lapply(labels, function(x) paste0(sort(x), collapse = "\n"))
#' 		cluster_xy <- split(data_plot[, 1:2], data_plot[, 3])
#' 		cluster_xy <- do.call(rbind, lapply(cluster_xy, colMeans))
#' 		for(i in 1:nrow(cluster_xy)){
#' 			text(cluster_xy[i, 1], cluster_xy[i, 2], labels[names(labels) == rownames(cluster_xy)[i]][[1]], cex=0.5, font=2)
#' 		}
#' 	}
#' 	
#' 	par(mar=c(0, 0, 0, 1))
#' 	plot.new()
#' 	legend("center", rev(col_levels), col=rev(pal), pch=16, bty="n", cex=0.6)
#' 	
#' 	
#' 	dev.off()
#' 	
#' 	
#' }
