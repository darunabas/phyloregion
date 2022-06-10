#' Conversion of community data
#'
#' These functions convert a community data to compressed sparse matrix,
#' dense matrix and long format (e.g. species records).
#'
#' @param x A community data which one wants to transform
#' @param grids column name of the column containing grid cells
#' @param species column name of the column containing  the species / taxa names
#' @rdname long2sparse
#' @importFrom Matrix sparseMatrix as.matrix
#'
#' @return A compressed sparse community matrix of sites by species
#'
#' @examples
#' data(africa)
#' africa$comm[1:5, 1:20]
#' long <- sparse2long(africa$comm)
#' long[1:5, ]
#' sparse <- long2sparse(long)
#' all.equal(africa$comm, sparse)
#'
#' dense_comm <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), paste0("sp", 1:4)))
#' dense_comm
#' sparse_comm <- dense2sparse(dense_comm)
#' sparse_comm
#' sparse2long(sparse_comm)
#'
#' @keywords Conversion
#' @export
long2sparse <- function(x, grids="grids", species="species"){
  x <- as.data.frame(x)
  x <- x[, c(grids, species)]
  names(x) <- c("grids", "species")
  grids <- factor(as.character(x$grids))
  species <- factor(as.character(x$species))
  res <- sparseMatrix(as.integer(grids), as.integer(species), x=1,
                      dimnames = list(levels(grids), levels(species)))
  res
}


#' @rdname long2sparse
#' @export
sparse2long <- function(x){
  if(!is(x, "sparseMatrix")) stop("I expected x to be a sparse matrix!")
  d <- as.data.frame(Matrix::summary(x))
  d$grids <- rownames(x)[d$i]
  d$species <- colnames(x)[d$j]
  if(!is.null(d$x)){
    ind <- rep(seq_along(d$x), d$x)
    dat <- d[ind, c("grids", "species")]
  } else dat <- d[, c("grids", "species")]
  dat
}


#' @rdname long2sparse
#' @export
dense2sparse <- function(x){
  x <- as.matrix(x)
  Matrix(x, sparse=TRUE)
}

#' @rdname long2sparse
#' @export
sparse2dense <- function(x){
  as.matrix(x)
}

#' @rdname long2sparse
#' @export
long2dense <- function(x){
   as.matrix(long2sparse(x))
}

#' @rdname long2sparse
#' @export
dense2long <- function(x){
  sparse2long(dense2sparse(x))
}

