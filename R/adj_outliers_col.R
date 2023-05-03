#' Remove outliers to avoid that these influence colors
#' @param x numeric vector
#' @description The definition of outliers is based on the function 'boxplot.stats()'. Values of x lower than the extreme of the lower whisker are set equal to it; values of x greater than the exteme of the upper whisker are set equal to it.
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
