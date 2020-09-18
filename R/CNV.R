#' cnv
#' @param x matrix for cnv calculation


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

