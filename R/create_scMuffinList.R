#' Create scMuffin list
#' @param dgCMatrix Accepted input: only dgCMatrix
#' @description Create scMuffin list
#' @export
#' 

create_scMuffinList <- function(dgCMatrix) {
    if (is(dgCMatrix, 'sparseMatrix') == TRUE) {
        scMuffinList$gene_by_cells <- as.data.frame(t(dgCMatrix))
    } else {
        warning("The only admitted input data is a dgCMatrix.")
    }
}


