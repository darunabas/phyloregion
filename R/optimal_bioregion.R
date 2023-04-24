.invalid <- function(x) {
  if (missing(x) || is.null(x) || length(x) == 0)
    return(TRUE)
  if (is.list(x))
    return(all(vapply(x, .invalid, FALSE)))
  else if (is.vector(x))
    return(all(is.na(x)))
  else return(FALSE)
}


.count.observation <- function(dist.obj) {
  nrow(as.matrix(dist.obj))
}


.PairwiseDistanceMatrixToVariance <- function(dist.obj) {
  mean(dist.obj**2) / 2
}


.retrieveDistanceMatrix <- function(dist.obj, indexVector) {
  as.dist(as.matrix(dist.obj)[indexVector, indexVector])
}


.tssByDistanceMatrix <- function(dist.obj) {
  n <- nrow(as.matrix(dist.obj))
  .PairwiseDistanceMatrixToVariance(dist.obj) * (n - 1)
}


.wssByDistanceMatrix <- function(dist.obj, indexVector) {
  n <- length(indexVector)
  if (n > 1) {
    i.dist.obj <- .retrieveDistanceMatrix(dist.obj, indexVector)
    .PairwiseDistanceMatrixToVariance(i.dist.obj) * (n - 1)
  } else if (n == 1) {
    return(0)
  } else {
    return(NA)
  }
}


.totwssByDistanceMatrix <- function(dist.obj, clusters) {
  cluster <- c()
  wss <- c()
  totwss <- 0
  for (i.cluster in unique(clusters)) {
    indexVector <- which(clusters == i.cluster)
    i.wss <- .wssByDistanceMatrix(dist.obj, indexVector)
    cluster <- c(cluster, i.cluster)
    wss <- c(wss, i.wss)
    ##     print(cluster)
    ##     print(indexVector)
    ##     print(i.wss)
  }
  totwss <- sum(wss)
  return(list(cluster = cluster, wss = wss, totwss = totwss))
}


.bssByDistanceMatrix <- function(dist.obj, clusters) {
  totwss <- .totwssByDistanceMatrix(dist.obj, clusters)$totwss
  tss <- .tssByDistanceMatrix(dist.obj)
  totbss <- tss - totwss
  ##   print(tss)
  ##   print(bss)
  return(list(totbss = totbss))
}

css <- function(dist.obj, clusters) {
  if (!(is.integer(clusters) & length(clusters) == .count.observation(dist.obj))) {
    stop("`clusters' should be a vector of integers of the same length as of the number of observations in `dist.obj'.")
  }
  ret <- list()
  ret$k <- length(unique(clusters))
  ret$wss <- .totwssByDistanceMatrix(dist.obj, clusters)$wss
  ret$totwss <- .totwssByDistanceMatrix(dist.obj, clusters)$totwss
  ret$totbss <- .bssByDistanceMatrix(dist.obj, clusters)$totbss
  ret$tss <- .tssByDistanceMatrix(dist.obj)
  meta <- list()
  meta$dist.obj <- dist.obj
  meta$clusters <- clusters
  attr(ret, "meta") <- meta
  class(ret) <- c("css", class(ret))
  return(ret)
}



.css.hclustS <- function(dist.obj, hclust.obj = NULL, hclust.FUN = hclust,
                         hclust.FUN.MoreArgs = list(method = "ward"), k) {
  #    if (.invalid(hclust.obj)) {
  #      hclust.obj <- .call.FUN(hclust.FUN,dist.obj,hclust.FUN.MoreArgs)
  #    }
  tmp.cutree <- cutree(hclust.obj, k = k)

  res <- css(dist.obj, clusters = tmp.cutree)
  attr(res, "meta")$hclust.obj <- hclust.obj
  res
}


.css.hclustM <- function(k, ...) {
  kVector <- k
  ret <- list()
  for (i in seq_along(kVector)) {
    ret[[i]] <- .css.hclustS(k = kVector[i], ...)
  }
  attr(ret, "meta")$hclust.obj <- attr(ret[[1]], "meta")$hclust.obj
  ret
}


css.hclust <- function(dist.obj,
                       hclust.obj = NULL,
                       hclust.FUN = hclust,
                       hclust.FUN.MoreArgs = list(method = "ward"),
                       k = NULL) {
  if (.invalid(k)) k <- min(.count.observation(dist.obj), 20)
  res <- .css.hclustM(dist.obj, hclust.obj, hclust.FUN, hclust.FUN.MoreArgs,
    k = 1:k)
  tss <- unlist(lapply(res, "[[", "tss"))
  totbss <- unlist(lapply(res, "[[", "totbss"))
  hclust.obj <-  attr(res, "meta")$hclust.obj
  ev <- totbss / tss
  ret <- data.frame(k = 1:k, ev = ev, totbss = totbss, tss = tss)

  meta <- list()
  meta$cmethod <- "hclust"
  meta$dist.obj <- dist.obj
  meta$k <- k
  meta$hclust.obj <- hclust.obj

  attr(ret, "meta") <- meta
  class(ret) <- c("css.multi", class(ret))

  ret
}

elbow <- function(x, inc.thres, ev.thres, precision = 3, print.warning = TRUE) {
  if (!inherits(x, "css.multi")) {
    stop("`x' should be an object of class `css.multi' from package GMD.")
  }

  if (inc.thres > 1 | inc.thres < 0) {
    stop("`inc.thres' should be a percentage value, ranging from 0 to 1.")
  }

  if (ev.thres > 1 | ev.thres < 0) {
    stop("`ev.thres' should be a percentage value, ranging from 0 to 1.")
  }

  if (length(inc.thres) > 1) {
    if (print.warning)
      warning("`inc.thres' has a length more than one; the 1st value is used.")
    inc.thres <- inc.thres[1]
  }

  if (length(ev.thres) > 1) {
    if (print.warning)
      warning("`ev.thres' has a length more than one; the 1st value is used.")
    ev.thres <- ev.thres[1]
  }

  ## increment (the first derivative) of totbss
  x$inc <- c(x$ev[-1] - x$ev[-nrow(x)], NA)
  ## x$dev1 <- x$inc
  ## x$dev2 <- c(x$dev1[-1]-x$dev1[-nrow(x)], NA)

  x$inc <- c(x$ev[-1] - x$ev[-nrow(x)], 0)

  ## print(x)
  ## cat(sprintf("x$ev=%s,x$inc=%s\n",x$ev,x$inc))
  k1 <- x[x$ev - ev.thres >= -0.1**precision, "k"][1]
  k2 <- x[-nrow(x), ][inc.thres - x[-nrow(x), ]$inc >= -0.1**precision, "k"][1]

  if (is.na(k1)) {
    if (print.warning)
      warning(sprintf("The explained variances is smaller than `ev.thres' when k=c(%s:%s); increase `k' or decrease `ev.thres' would help.", 1, max(x$k)))
  }
  if (is.na(k2)) {
    if (print.warning)
      warning(sprintf("The increment of explained variance is larger than `inc.thres' when k=c(%s:%s); increase `k' or increase `inc.thres' would help.", 1, max(x$k)))
  }

  k <- max(k1, k2)
  ret <- list(k = k, ev = x$ev[x$k == k], inc.thres = inc.thres,
              ev.thres = ev.thres)
  class(ret) <- c("elbow", class(ret))
  return(ret)
}



elbow.batch <- function(x, inc.thres = c(0.01, 0.05, 0.1),
                        ev.thres = c(0.95, 0.9, 0.8, 0.75, 0.67, 0.5, 0.33),
                        precision = 3) {
  ret <- NULL
  inc.thres <- sort(inc.thres, decreasing = FALSE)
  ev.thres <- sort(ev.thres, decreasing = TRUE)
  for (i.ev.thres in ev.thres) {
    for (i.inc.thres in inc.thres) {
      res <- elbow(x, inc.thres = i.inc.thres, ev.thres = i.ev.thres,
                   precision = precision, print.warning = FALSE)
      k <- res$k
      if (!is.na(k)) {
        info <- sprintf("A \"good\" k=%s (EV=%.2f) is detected when the EV is no less than %s\nand the increment of EV is no more than %s for a bigger k.\n", k, res$ev, i.ev.thres, i.inc.thres)
        ret <- list(k = k, ev = res$ev, inc.thres = i.inc.thres,
                    ev.thres = i.ev.thres)
        attr(ret, "description") <- info
        class(ret) <- c("elbow", class(ret))
        return(ret)
      }
    }
  }
  if (is.null(ret)) {
    stop("A good `k' is not available with provided inc.thres and ev.thres; please make adjustment, e.g. decrease `ev.thres', increase `inc.thres' or increase `k'.")
  }
}


#' Determine optimal number of clusters
#'
#' This function divides the hierarchical dendrogram into meaningful
#' clusters ("phyloregions"), based on the ‘elbow’ or ‘knee’ of
#' an evaluation graph that corresponds to the point of optimal curvature.
#
#' @param x a numeric matrix, data frame or \dQuote{dist} object.
#' @param method the agglomeration method to be used. This should
#' be (an unambiguous abbreviation of) one of \dQuote{ward.D}, \dQuote{ward.D2},
#' \dQuote{single}, \dQuote{complete}, \dQuote{average} (= UPGMA),
#' \dQuote{mcquitty} (= WPGMA),
#' \dQuote{median} (= WPGMC) or \dQuote{centroid} (= UPGMC).
#' @param k numeric, the upper bound of the number of clusters to
#' compute. DEFAULT: 20 or the number of observations (if less than 20).
#' @keywords phyloregion
#' @importFrom stats hclust as.dist
#'
#' @references
#' Salvador, S. & Chan, P. (2004) \emph{Determining the number of
#' clusters/segments in hierarchical clustering/segmentation algorithms}.
#' Proceedings of the Sixteenth IEEE International Conference on Tools
#' with Artificial Intelligence, pp. 576–584. Institute of Electrical
#' and Electronics Engineers, Piscataway, New Jersey, USA.
#'
#' Zhao, X., Valen, E., Parker, B.J. & Sandelin, A. (2011) Systematic
#' clustering of transcription start site landscapes.
#' \emph{PLoS ONE} \strong{6}: e23409.
#'
#' @return a list containing the following as returned from the GMD
#' package (Zhao et al. 2011):
#' \itemize{
#'   \item \code{k}:	optimal number of clusters (bioregions)
#'   \item \code{totbss}:	total between-cluster sum-of-square
#'   \item \code{tss}:	total sum of squares of the data
#'   \item \code{ev}:	explained variance given k
#' }
#' @examples
#' data(africa)
#' tree <- africa$phylo
#' bc <- beta_diss(africa$comm)
#' (d <- optimal_phyloregion(bc[[1]], k=15))
#' plot(d$df$k, d$df$ev, ylab = "Explained variances",
#'   xlab = "Number of clusters")
#' lines(d$df$k[order(d$df$k)], d$df$ev[order(d$df$k)], pch = 1)
#' points(d$optimal$k, d$optimal$ev, pch = 21, bg = "red", cex = 3)
#' points(d$optimal$k, d$optimal$ev, pch = 21, bg = "red", type = "h")
#' @export
optimal_phyloregion <- function(x, method = "average", k = 20) {
  m <- hclust(as.dist(x), method = method)
  m <- css.hclust(as.dist(x), m, k = k)
  n <- elbow.batch(m)
  return(list(df = m, optimal = n))
}
