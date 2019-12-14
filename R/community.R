#' Read in sparse community matrices
#'
#' read.community reads in file containing occurence data and returns a sparse
#' matrix.
#'
#' @param file A file name.
#' @param ... further arguments passed to or from other methods.
#' @keywords cluster
#' @rdname read.community
#' @importFrom ape Nnode Ntip
#' @importFrom Matrix Matrix sparseMatrix t
#' @importFrom phangorn as.splits
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


as.splits.matrix <- function(x){
  y <- vector("list", nrow(x))
  tmp <- seq_len(ncol(x))
  for(i in seq_len(nrow(x))){
    y[[i]] <- tmp[x[i,]>0]
  }
  names(y) <- rownames(x)
  attr(y, "labels") <- colnames(x)
  class(y) <- c("splits")
  y
}

as.splits.data.frame <- function(x){
  as.splits(as.matrix(x))
}


as.splits.Matrix <- function(x, ...){
#  if(is.data.frame(x)) x <- as.matrix(x)
#  if(is.matrix(x)) x <- Matrix(x, sparse=TRUE, doDiag = FALSE)
  M <- t(x)
  res <- vector("list", ncol(M))
  for(i in seq_len(ncol(M))){
    res[[i]] <- M@i[(M@p[i]+1) : (M@p[i+1])] + 1
  }
  names(res) <- colnames(M)
  attr(res, "labels") <- row.names(M)
  class(res) <- c("splits")
  res
}


