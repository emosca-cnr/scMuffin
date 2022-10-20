#' Over Representation Analysis
#' @param wb hits (white balls)
#' @param bb other elements (black balls)
#' @param gsl named list of sets
#' @param p_adj_method p value adjustment method, see p.adjust.methods
#' @importFrom stats p.adjust
#' @description Over representation analysis
#' @export

ora <- function(wb, bb, gsl, p_adj_method='fdr'){


  out <- lapply(gsl, function(x) ora1gs(wb, bb, x))

  out <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
  out$N <- length(wb) + length(bb)
  out$exp <- out$wb * out$bd / out$N
  out$id <- rownames(out)
  out$p_adj <- stats::p.adjust(out$p, method = p_adj_method)
  out$er <- out$wbd / out$exp

  return(out[, c('id', 'N', 'wb', 'bb', 'bd', 'wbd', 'exp', 'er', 'p', 'p_adj')])

}
