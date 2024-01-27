JENKS <- function(x, k){
  d <- sort(x)
  k <- k
  mat1 <- matrix(1, length(d), k)
  mat2 <- matrix(0, length(d), k)
  mat2[2:length(d), 1:k] <- .Machine$double.xmax
  v <- 0
  for (l in 2:length(d)) {
    s1 = s2 = w = 0
    for (m in 1:l) {
      i3 <- l - m + 1
      val <- d[i3]
      s2 <- s2 + val * val
      s1 <- s1 + val
      w <- w + 1
      v <- s2 - (s1 * s1)/w
      i4 <- trunc(i3 - 1)
      if (i4 != 0) {
        for (j in 2:k) {
          if (mat2[l, j] >= (v + mat2[i4, j - 1])) {
            mat1[l, j] <- i3
            mat2[l, j] <- v + mat2[i4, j - 1]
          }
        }
      }
    }
    mat1[l, 1] <- 1
    mat2[l, 1] <- v
  }
  kclass <- 1:k
  kclass[k] <- length(d)
  k <- length(d)
  last <- length(d)
  for (j in length(kclass):1) {
    id <- trunc(mat1[k, j]) - 1
    kclass[j - 1] <- id
    k <- id
    last <- k - 1
  }
  brks <- d[c(1, kclass)]
}

#' Bin values
#'
#' \code{choropleth} discretizes the values of a quantity for mapping.
#'
#' @param x Vector of values to discretize.
#' @param k Numeric, the desired number of bins to discretize.
#' @param breaks one of \dQuote{equal}, \dQuote{pretty}, \dQuote{jenks},
#' \dQuote{quantile} or numeric vector with the actual breaks by
#' specifying the minimum (\code{min}) and maximum (\code{max}) bounds.
#' @param min the minima of the lowest bound of the break.
#' @param max the maxima of the upper bound of the break
#' @rdname choropleth
#' @keywords bioregion
#' @importFrom stats quantile
#' @importFrom terra classify values crop mask levels
#' @return a vector with the discretized values.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#' @seealso \code{\link{coldspots}}
#' @export
choropleth <- function(x, k = 10, breaks = "quantile", min = NULL, max = NULL) {
  ras <- x
  nm <- names(ras)
  x[x==0] <- NA
  r <- values(x)
  
  if (!is.null(min) & !is.null(max)) {
    y <- seq(min, max, length.out = length(r))
    quants <- switch(breaks,
                     quantile = quantile(y, probs = seq(0, 1, 1/k)),
                     equal = seq(min(y), max(y), length.out = (k + 1)),
                     pretty = c(pretty(y, k = k)),
                     jenks = JENKS(y, k))
  } else
    quants <- switch(breaks,
                     quantile = quantile(r |> na.omit(), probs=seq(0, 1, 1/k)),
                     equal = seq(min(r |> na.omit()), max(r |> na.omit()), 
                                 length.out = (k + 1)),
                     pretty = c(pretty(r |> na.omit(), k = k)),
                     jenks = JENKS(r |> na.omit(), k))
  if(length(which(quants==0)) > 0) {
    brks <- unique(c(0, ceiling(quants[-which(quants == 0)])))
  }else{
    brks <- unique(c(0, ceiling(quants)))
  }
  y <- classify(x, c(0, round(brks)))
  levels(y)[[1]][[nm]] <- brks
  y[is.na(y)] <- 0
  z <- y |> crop(ras) |> mask(ras)
  return(z)
}

