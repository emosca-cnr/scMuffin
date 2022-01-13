#' gene_set_score
#' @param gene_set gene set
#' @param genes_by_cells real data
#' @param bins random data
#' @param nmark_min min number of markers
#' @param ncells_min min number of cells
#' @param k number of random gene set to be considered
#' @param kmin minimim number of k
#' @param verbose verbosity
#' @param null_model boolean, whether to use the permutations or not
#' @return list with
#'    score_table: data frame with gene set score for each cell
#'    permutations: list with the permutations
#' @author Ettore Mosca
#' @import graphics
#' @export

#gene_set_score <- function(gene_set_data, control_set_data, nmark_min=5, ncells_min=NULL){
gene_set_score <- function(gene_set, genes_by_cells, bins, nmark_min=5, ncells_min=5, k=100, kmin=50, verbose=TRUE, null_model=TRUE){
	
	if(sum(rownames(genes_by_cells) %in% gene_set) < nmark_min){
		message("not enough genes in gene set, please consider only gene set with enough number of genes\n")
		res <- data.frame(case=rep(NA, ncol(genes_by_cells)), case.N=NA, case.AV=NA, nmark_min=F, avg_control=NA, control.AV=NA, null_ok=F, avg_delta_score=NA, delta_score=NA, stringsAsFactors = F, row.names = colnames(genes_by_cells))
	}else{
		gene_set_data <- genes_by_cells[rownames(genes_by_cells) %in% gene_set, , drop=F]
		gene_set_data[gene_set_data==0] <- NA #NA to avoid inclusion in mean calculation
		
		if(null_model){
			control_set_data <- sc_create_null(genes_by_cells, bins = bins, gene_set = gene_set, k = k)
			cat("nulls obtained: ", length(control_set_data), "\n")
			for(i in 1:length(control_set_data)){
				control_set_data[[i]][control_set_data[[i]]==0] <- NA
			}
			
			#null_mod_gs <- list(gs=hist(log2(rowSums(gene_set_data, na.rm = T)), 20, plot = F))
			#null_mod_gs <- c(null_mod_gs, lapply(control_set_data, function(x) hist(log2(rowSums(x, na.rm = T)), 20, plot = F)))
			
		}
		
		#old version
		#ans <- data.frame(case=colMeans(gene_set_data), control=colMeans(control_set_data), case.N=nrow(gene_set_data), case.AV=colSums(gene_set_data!=0), case.NA=colSums(gene_set_data==0), control.N=nrow(control_set_data), control.AV=colSums(control_set_data != 0), control.NA=colSums(control_set_data==0), stringsAsFactors = F)
		#ans$score <- ans$case - ans$control
		
		#new version
		ans <- list(gs=data.frame(case=colMeans(gene_set_data, na.rm = T), case.N=nrow(gene_set_data), case.AV=colSums(!is.na(gene_set_data)), stringsAsFactors = F)) #real
		
		#z-scores....
		#temp <- gene_set_data
		#for(i in 1:nrow(temp)){
		#  temp_idx <- !is.na(temp[i, ])
		#  temp[i, temp_idx] <- scale(temp[i, temp_idx])
		#}
		
		#ans[[1]] <- cbind(ans[[1]], case.z=colMeans(temp, na.rm = T))
		
		if(null_model){
			ans <- c(ans, lapply(control_set_data, function(x) data.frame(control=colMeans(x, na.rm = T), control.N=nrow(x), control.AV=colSums(!is.na(x)), stringsAsFactors = F))) #controls
		}
		
		#case and control will have NaN for cells with all NA
		for(i in 1:length(ans)){
			ans[[i]][is.nan(ans[[i]][, 1]), 1] <- NA
		}
		
		ans <- lapply(ans, function(x) cbind(x, nmark_min=x[,3]>=nmark_min)) #only cells in which the gene set is available
		
		if(null_model){
			#cells with enough genes
			vect_selected <- do.call(cbind, t(lapply(ans, function(x) x$nmark_min))) #the first column contains real data
			
			#null_ok <- rowSums(vect_selected[, -1])  > (ncol(vect_selected[, -1])) #old
			gene_set_ok_in_nulls <- rowSums(vect_selected[, -1], na.rm = T) #number of permutations in which each cell satistfies nmark_min
			if(verbose){
				cat("gene_set_ok_in_nulls...\n")
				print(table(gene_set_ok_in_nulls[gene_set_ok_in_nulls>0]))
			}
			
			if(max(gene_set_ok_in_nulls) < kmin){
				cat("kmin", kmin, "; found, ", max(gene_set_ok_in_nulls), "\n")
				cat("Consider reducing kmin\n")
			}
			#null_ok <- gene_set_ok_in_nulls > (ncol(vect_selected)-1)*0.9 #nubmer of cells that have a sufficient number of available permutations
			null_ok <- gene_set_ok_in_nulls >= kmin
			
			
			if(sum(vect_selected[, 1]) >= ncells_min){
				
				vect_selected[, 1] <- vect_selected[, 1] & null_ok
				
				if(sum(vect_selected[, 1]) >= ncells_min){
					
					#score matrix
					score <- do.call(cbind, t(lapply(ans[-1], function(x) ans[[1]]$case - x$control)))
					#score <- score * sign(vect_selected[, -1]) #remove permutations not valid
					#score[score==0] <- NA #0 appeared due sign
					score[!vect_selected[, -1]] <- NA #new
					score <- rowMeans(score, na.rm = T)
					score[is.nan(score)] <- NA
					
					#average control value for each cell
					avg_control <- do.call(cbind, t(lapply(ans[-1], function(x) x$control)))
					#avg_control <- avg_control * sign(vect_selected[, -1])
					#avg_control[avg_control==0] <- NA
					avg_control[!vect_selected[, -1]] <- NA #new
					avg_control <- rowMeans(avg_control, na.rm = T)
					avg_control[is.nan(avg_control)] <- NA
					
					res <- data.frame(ans[[1]], avg_control=avg_control, control.AV=gene_set_ok_in_nulls, null_ok=null_ok, avg_delta_score=score, delta_score=ans[[1]]$case-avg_control, stringsAsFactors = F) #difference between real gs score (1), and each of the null values [2, k]
					
					#res <- res[res$nmark_min, ] #eliminate cells without a sufficient number of genes
					
				}else{
					message("cannot create null models for at least ", ncells_min, " cells\n")
					res <- data.frame(ans[[1]], avg_control=NA, control.AV=NA, null_ok=F, avg_delta_score=NA, delta_score=NA, stringsAsFactors = F)
				}
				
				#res <- res[res$nmark_min, ] #eliminate cells without a sufficient number of genes
				
			}else{
				
				message("cannot calcolate gene set score for at least ", ncells_min, " cells\n")
				res <- data.frame(case=rep(NA, ncol(genes_by_cells)), case.N=NA, case.AV=NA, nmark_min=F, avg_control=NA, control.AV=NA, null_ok=F, avg_delta_score=NA, delta_score=NA, stringsAsFactors = F, row.names = colnames(genes_by_cells))
				
			}
			
		}else{ #without nullmodel
			
			res <- data.frame(ans[[1]], avg_control=NA, control.AV=NA, null_ok=NA, avg_delta_score=NA, delta_score=NA, stringsAsFactors = F)
			
		}
		#ans <- list(score_table=res, permutations=ans[-1])
	}
	return(res)
	
}
