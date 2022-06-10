## This file is borrowed from Kurt Hornik's relations package

## A simple class for sparse (triplet) matrices.

## Mostly intended for being able to take advantage of LP solvers which
## allow for sparse specifictions of (possible rather large) constraint
## matrices.

simple_triplet_matrix <-
    function(i, j, v, nrow = max(i), ncol = max(j))
    {
        structure(list(i = i, j = j, v = v, nrow = nrow, ncol = ncol),
                  class = "simple_triplet_matrix")
    }

as.simple_triplet_matrix <-
    function(x)
        UseMethod("as.simple_triplet_matrix")

as.simple_triplet_matrix.simple_triplet_matrix <- identity

as.simple_triplet_matrix.matrix <-
    function(x)
    {
        if(prod(dim(x)) == 0L)
            return(simple_triplet_matrix(integer(), integer(), c(x),
                                         nrow = nrow(x), ncol = ncol(x)))
        ind <- which(x != vector(typeof(x), 1L), arr.ind = TRUE)
        simple_triplet_matrix(ind[, 1L], ind[, 2L], x[ind],
                              nrow = nrow(x), ncol = ncol(x))
    }

as.matrix.simple_triplet_matrix <-
    function(x, ...)
    {
        nr <- x$nrow
        nc <- x$ncol
        y <- matrix(vector(typeof(x$v), nr * nc), nr, nc)
        y[cbind(x$i, x$j)] <- x$v
        y
    }

## We could also simply write a method to coerce to a dgTMatrix, based
## on something like
##  new("dgTMatrix",
##       i = as.integer(i - 1),
##       j = as.integer(j - 1),
##       x = v,
##       Dim = c(nrow, ncol))
## (Note that these have C-style indices starting at zero.)

dim.simple_triplet_matrix <-
    function(x)
        c(x$nrow, x$ncol)

`[.simple_triplet_matrix` <-
    function(x, i, j, ...)
    {
        ## (Well, we certainly don't drop ...)
        mi <- missing(i)
        mj <- missing(j)
        na <- nargs()

        if(mi && mj) {
            out <- vector(mode = typeof(x$v), length = x$nrow * x$ncol)
            out[(x$j - 1) * x$nrow + x$i] <- x$v
        }
        else if(na == 2L) {
            ## Single index subscripting.
            if(is.logical(i))
                stop("Logical subscripting currently not implemented.")
            else if(is.character(i))
                stop("Character subscripting currently not implemented.")
            else if(!is.matrix(i)) {
                ## Let's hope we have a vector.
                ## What if we have both negatives and positives?
                if(all(i >= 0)) {
                    i <- i[i > 0]
                    out <- vector(mode = typeof(x$v), length = length(i))
                    pos <- match(i, (x$j - 1) * x$nrow + x$i, 0)
                    out[i[pos > 0]] <- x$v[pos]
                    out
                } else if(all(i <= 0)) {
                    out <- vector(mode = typeof(x$v),
                                  length = x$nrow * x$ncol)
                    out[(x$j - 1) * x$nrow + x$i] <- x$v
                    out <- out[i]
                }
                else stop("Cannot mix positive and negative subscripts.")
            }
            else {
                ## Note that negative values are not allowed in a matrix
                ## subscript.
                if((ncol(i) != 2L) || (any(i < 0)))
                    stop("Invalid subscript.")
                i <- i[!apply(i == 0, 1L, any), , drop = FALSE]
                out <- vector(mode = typeof(x$v), length = nrow(i))
                pi <- match(i[, 1L], x$i)
                pj <- match(i[, 2L], x$j)
                ind <- which(pi == pj)
                out[ind] <- x$v[pi[ind]]
            }
        }
        else {
            nr <- x$nrow
            nc <- x$ncol
            ## Two index subscripting is rather tricky, as it can also be
            ## used for rearranging and "recycling" rows and columns.  Let
            ## us not support the latter for now, so that selected rows and
            ## columns must be unique.
            if(mi) {
                pos <- rep.int(TRUE, length(x$v))
                nr <- x$nrow
                pi <- seq_len(nr)
            }
            else if(!is.numeric(i))
                stop("Only numeric two-index subscripting is implemented.")
            else {
                pi <- seq_len(x$nrow)
                if(all(i >= 0)) {
                    i <- i[i > 0]
                    if(any(duplicated(i)))
                        stop("Repeated indices currently not allowed.")
                } else if(all(i <= 0))
                    i <- pi[i]
                else
                    stop("Cannot mix positive and negative subscripts.")
                nr <- length(i)
                pos <- match(x$i, i, 0) > 0
                pi[i] <- seq_len(nr)
            }
            if(mj) {
                nc <- x$ncol
                pj <- seq_len(nc)
            }
            else if(!is.numeric(j))
                stop("Only numeric 2-index subscripting is implemented.")
            else {
                pj <- seq_len(x$ncol)
                if(all(j >= 0)) {
                    j <- j[j > 0]
                    if(any(duplicated(j)))
                        stop("Repeated indices currently not allowed.")
                } else if(all(j <= 0))
                    j <- pj[j]
                else
                    stop("Cannot mix positive and negative subscripts.")
                nc <- length(j)
                pos <- pos & (match(x$j, j, 0) > 0)
                pj[j] <- seq_len(nc)
            }
            out <- simple_triplet_matrix(pi[x$i[pos]],
                                         pj[x$j[pos]],
                                         x$v[pos], nr, nc)
        }

        out
    }

rbind.simple_triplet_matrix <-
    function(..., deparse.level = 1L)
    {
        ## Ignore 'deparse.level' ...
        Reduce(function(x, y) {
            if((nc <- ncol(x)) != ncol(y))
                stop("Numbers of columns of matrices must match.")
            nr <- nrow(x)
            simple_triplet_matrix(c(x$i, y$i + nr),
                                  c(x$j, y$j),
                                  c(x$v, y$v),
                                  nrow = nr + nrow(y), ncol = nc)
        },
        lapply(Filter(Negate(is.null), list(...)),
               as.simple_triplet_matrix))
    }

cbind.simple_triplet_matrix <-
    function(..., deparse.level = 1L)
    {
        ## Ignore 'deparse.level' ...
        Reduce(function(x, y) {
            if((nr <- nrow(x)) != nrow(y))
                stop("Numbers of rows of matrices must match.")
            nc <- ncol(x)
            simple_triplet_matrix(c(x$i, y$i),
                                  c(x$j, y$j + nc),
                                  c(x$v, y$v),
                                  nrow = nr, ncol = nc + ncol(y))
        },
        lapply(Filter(Negate(is.null), list(...)),
               as.simple_triplet_matrix))
    }

t.simple_triplet_matrix <-
    function(x)
        simple_triplet_matrix(x$j, x$i, x$v, x$ncol, x$nrow)

## Utitilies for creating special simple triplet matrices:

.simple_triplet_zero_matrix <-
    function(nrow, ncol = nrow, mode = "double")
        simple_triplet_matrix(integer(), integer(), vector(mode, 0L),
                              nrow, ncol)

.simple_triplet_diag_matrix <-
    function(x, nrow)
    {
        x <- rep(x, length.out = nrow)
        i <- seq_len(nrow)
        simple_triplet_matrix(i, i, x, nrow, nrow)
    }

#' @importFrom Matrix Matrix
#' @export
Matrix::Matrix

#' @importFrom Matrix sparseMatrix
#' @export
Matrix::sparseMatrix

#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
