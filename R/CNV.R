#' CNV
#' @param x list of dataframes retrieved by 'preprocess_object_for_cnv'.
#' @description Function to compute CNV.
#' @details As a list of dataframes, it is convenient to use 'mclapply' from parallel package. 
#' @usage mclapply(mat_sorted, function(x) apply(x[,5:dim(x)[2]], 2, function( y ) CNV(y)), mc.cores = 2)
#' @references Patel paper supplementary material (?)
#' @author Valentina Nale


CNV <- function(x) {
  
  N <- length(x)

  Ek <- c()
  
  # first check: length
  if (N<101) {
    stop('The number of gene is not sufficient to compute CNV')
  } else {
    symbol_names <- names(x)[51:(N-51)]
    for (i in 51:(N-51)) {
      val <- mean(x[(i-50):(i+50)])
      Ek <- c(Ek, val)
      Ek <<- Ek
      }
  }
  names(Ek)<-symbol_names
  return(Ek)
}

