#' Convert community matrix to sparse community matrix
#' 
#' Converts a species by sites community matrix to a saprse community matrix. 
#' Expects a matrix with rows as sites and columns as species or taxa.
#' 
#' @param x a community matrix or data frame
#' 
#' @note See \url{https://github.com/darunabas/bioregion} for more
#' details and tutorials.
#' 
#' @importFrom methods as
#' @keywords Conversion
#' @examples
#' require(sp)
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), paste0("sp", 1:4)))
#' comm2sparse(com)
#' @export
comm2sparse <- function(x){
  x <- as.matrix(x)
  m <- as(x, "ngCMatrix")
  m
}

