#' Mean distance matrix from a set of distance matrices
#'
#' This function generates the mean pairwise distance matrix from a set
#' many pairwise distance matrices. Note: all matrices should be of the same dimension.
#'
#' @param files list of pairwise distance matrices stored as CSVs or .rds
#' with the same dimensions.
#' @param trace Trace the function; trace = 2 or higher will be more voluminous.
#' @param ... Further arguments passed to or from other methods.
#' @rdname mean_dist
#' @return average distance matrix
#' @importFrom utils read.csv txtProgressBar setTxtProgressBar
#'
#' @export
mean_dist <- function(files, trace = 1, ...){
  tips <- sort(labels(files[[1]]))
  ntips <- length(tips)
  res <- matrix(0, ntips, ntips, dimnames = list(tips, tips))
  tmp <- matrix(0L, ntips, ntips, dimnames = list(tips, tips))
  if (interactive() && trace > 0) {
    pb <- txtProgressBar(min = 0, max = length(files), style = 3,
                         width = getOption("width")/2L)
  }
  for(i in seq_along(files)){
    d <- as.matrix(files[[i]])
    dnam <- colnames(d)
    res[dnam,dnam] <- res[dnam,dnam] + d
    tmp[dnam,dnam] <- tmp[dnam,dnam] + 1L
    if (interactive() && trace > 0) setTxtProgressBar(pb, i)
  }
  res[tmp>0L] <- res[tmp>0L] / tmp[tmp>0L]
  res <- as.dist(res)
  res
}

