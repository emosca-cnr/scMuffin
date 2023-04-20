#' Remove outliers to avoid that these influence colors
#' @param x numeric vector
#' @description Remove outliers to avoid that these influence colors
#' @export

adj_outliers_col <- function(x){
  
	#print(summary(x))
  s <- boxplot.stats(x)$stats
  #print(s)
  
  ans <- x
  ans[ans < s[1] ] <- s[1]
  ans[ans > s[5] ] <- s[5]
  #print(summary(ans))
  
  return(ans)
  
}
