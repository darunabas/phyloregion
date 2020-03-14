#' Mean distance matrix from a set of distance matrices
#'
#' This function generates the mean pairwise distance matrix from a set
#' many pairwise distance matrices. Note: all matrices should be of the same dimension.
#'
#' @param files list of pairwise distance matrices stored as CSVs or .rds
#' with the same dimensions.
#' @param tips list of site or grid names
#' @param ... Further arguments passed to or from other methods.
#' @rdname mean_dist
#' @return average distance matrix
#' @importFrom utils read.csv
#'
#' @export
mean_dist <- function(files, tips, ...){
  ntips <- length(tips)
  res <- matrix(0, ntips, ntips, dimnames = list(tips, tips))
  tmp <- matrix(0L, ntips, ntips, dimnames = list(tips, tips))
  for(i in seq_along(files)){
    #d <-  read.csv(files[i], row.names = 1)
    d <- as.matrix(files[[i]])
    dnam <- colnames(d)
    res[dnam,dnam] <- res[dnam,dnam] + d
    tmp[dnam,dnam] <- tmp[dnam,dnam] + 1L
  }
  res[tmp>0L] <- res[tmp>0L] / tmp[tmp>0L]
  res <- as.dist(res)
  res
}

