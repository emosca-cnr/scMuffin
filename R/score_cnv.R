#' Score CNV
#' @param cnv_res Object resulting from calculate_cnv function
#' @param plot_score NULL. Set as TRUE it return a plot.
#' @export

score_cnv <- function(cnv_res, plot_score=FALSE) {
  
	X <- cnv_res
	
  #### 1) sum of squares of CNV
  cs <- colSums(X^2)

  ##### 2) correlation of each cell CNV vector with the average CNV vector of top 10% of cells

  top10 <- quantile(cs, probs=0.9)
  cs_filt <- names(cs[cs>top10])
  cnv_avg_top10 <- rowMeans(X[, colnames(X) %in% cs_filt])


  corr <- cor(cnv_res, cnv_avg_top10, method = "spearman")
  
  res_corr <- as.data.frame(corr)
  res_cs <- as.data.frame(cs)
  res <- merge(res_cs, res_corr, by=0, all=TRUE)
  rownames(res) <- res$Row.names
  res$Row.names <- NULL
  names(res)[names(res) == "V1"] <- "CNV_R"
  names(res)[names(res) == "cs"] <- "CNV_signal"
  plot_score=TRUE

  if (plot_score) {
    grDevices::jpeg("plot_score.jpg", width = 180, height = 180, res = 300, units = "mm")
    plot(res$CNV_signal, res$CNV_R, xlab="CNV signal", ylab="CNV correlation")
    dev.off()
  }

  return(res)
}
