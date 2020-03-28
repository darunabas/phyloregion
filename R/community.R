#' Read in sparse community matrices
#'
#' \code{read.community} reads in file containing occurrence data and returns a
#' sparse matrix.
#'
#' @param file A file name.
#' @param grids Column name of the column containing grid cells.
#' @param species Column name of the column containing the species / taxa names.
#' @param ... further arguments passed to or from other methods.
#' @return \code{read.community} returns a sparse matrix (an object of class
#' "dgCMatrix").
#' @keywords cluster
#' @rdname read.community
#' @importFrom ape Nnode Ntip
#' @importFrom Matrix Matrix sparseMatrix t
#' @importFrom methods is
#' @importFrom utils read.csv
#' @examples
#' \donttest{
#' df <- data.frame(grids=paste0("g", c(1,1,2,3,3)),
#'                  species = paste0("sp", c(1,3,2,1,4)))
#' df
#' tmp <- tempfile()
#' write.csv(df, tmp)
#' (M <- read.community(tmp) )
#' sparse2long(M)
#' unlink(tmp)
#' }
#' @export
read.community <- function(file, grids="grids", species="species", ...){
  d <- read.csv(file, ...)
  M <- Matrix::sparseMatrix(as.integer(d[,grids]), as.integer(d[,species]),
                  x = rep(1L, nrow(d)),
                  dimnames = list(levels(d[,"grids"]), levels(d[,"species"])))
  M
}



