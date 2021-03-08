#' Calculate the overlap coefficient between two sets
#' @export
overlap_matrix <- function(x){
	
	all_states <- split(x[[1]], x[[1]])
	for(i in 2:length(x)){
		all_states <- c(all_states, split(x[[i]], x[[i]]))
	}
	all_states_n <- length(all_states)
	for(i in 1:all_states_n){
		all_states[[i]] <- names(all_states[[i]])
	}
	
	ans <- matrix(0, all_states_n, all_states_n, dimnames = list(names(all_states), names(all_states)))
	for(i in 2:all_states_n){
		for(j in 1:(i-1)){
			ans[i,j] <- overlap_coefficient(all_states[[i]], all_states[[j]])
		}
	}
	ans <- ans + t(ans)
	diag(ans) <- 1
	
	return(ans)
	
}