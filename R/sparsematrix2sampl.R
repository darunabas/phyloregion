#' Deconstruct sparse matrix to community data
#'
#' This function converts a compressed sparse matrix to a base R
#' community data frame of (row, column) pairs.
#'
#' @param x A compressed sparse matrix of class \dQuote{dgCMatrix}, or \dQuote{dgRMatrix}.
#' @rdname sparse2sampl
#' @keywords phyloregion
#' @importFrom Matrix summary
#'
#' @return a community data frame of sites and species pairs
#'
#' @examples
#' data(africa)
#' tree <- africa$phylo
#' M <- sampl2sparse(africa$comm, method = "nonphylo")
#' N <- sparse2sampl(M)
#' @export sparse2sampl
sparse2sampl <- function(x){
  d <- as.data.frame(Matrix::summary(x))
  d$grids <- rownames(x)[d$i]
  d$species <- colnames(x)[d$j]
  if(!is.null(d$x)){
    ind <- rep(seq_along(d$x), d$x)
    dat <- d[ind, c("grids", "species")]
  } else dat <- d[, c("grids", "species")]
  dat
}




