#' Mean matrix from a set of pairwise distance matrices
#'
#' This function generates the mean pairwise distance matrix from a set
#' many pairwise distance matrices. Note: all matrices should be of the same dimension.
#'
#' @param files list of pairwise distance matrices stored as CSVs or .rds
#' with the same dimensions.
#' @param tips list of site or grid names
#' @param verbose logical; if TRUE, show even more when running example code.
#' @param ... Further arguments passed to or from other methods.
#' @rdname av.matrix
#' @return average pairwise distance matrix
#' @importFrom utils read.csv
#'
#' @examples
#' rm(list = ls())
#' wd <- tempdir()
#' setwd(wd)
#'
#' ## Generate the sets of pairwise distance matrices
#' require(ape)
#' set.seed(1)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' tree <- rmtree(10, 4)
#' com <- matrix(c(1,0,1,1,0,0,
#'                1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), paste0("t", 1:4)))
#'
#'
#' for (i in 1:length(tree)) {
#'   pbc <- phylobeta(com, tree[[i]])
#'   write.csv(pbc, file = xzfile(paste0("matrix", i, ".xz")))
#' }
#'
## compute average pairwise distance matrix
#' files <- list.files(wd, pattern="*.xz", full.names=TRUE)
#' tr <- colnames(read.csv(files[1], row.names = 1, stringsAsFactors = FALSE))
#' tr1 <- sort(tr)
#' av.matrix(files, tr1)
#' @export
av.matrix <- function(files, tips, verbose = TRUE, ...){
  if (verbose) {
    message("Generating average pairwise matrix from several matrices")
  }
  ntips <- length(tips)
  res <- matrix(0, ntips, ntips, dimnames = list(tips, tips))
  tmp <- matrix(0L, ntips, ntips, dimnames = list(tips, tips))
  for(i in seq_along(files)){
    d <-  read.csv(files[i], row.names = 1)
    d <- as.matrix(d)
    dnam <- colnames(d)
    res[dnam,dnam] <- res[dnam,dnam] + d
    tmp[dnam,dnam] <- tmp[dnam,dnam] + 1L
  }
  res[tmp>0L] <- res[tmp>0L] / tmp[tmp>0L]
  res
}

