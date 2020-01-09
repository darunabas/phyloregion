#' Read in sparse community matrices
#'
#' \code{read.community} reads in file containing occurence data and returns a
#' sparse matrix.
#'
#' @param file A file name.
#' @param ... further arguments passed to or from other methods.
#' @keywords cluster
#' @rdname read.community
#' @importFrom ape Nnode Ntip
#' @importFrom Matrix Matrix sparseMatrix t
#' @importFrom methods is
#' @importFrom utils read.csv
#' @examples
#' data(africa)
#' @export
read.community <- function(file, ...){
  d <- read.csv(file, ...)
  M <- Matrix::sparseMatrix(as.integer(d[,"grids"]), as.integer(d[,"species"]),
                  x = rep(1L, nrow(d)),
                  dimnames = list(levels(d[,"grids"]), levels(d[,"species"])))
  M
}



