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
#' @param style one of \dQuote{equal}, \dQuote{pretty}, or \dQuote{quantile}.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname choropleth
#' @keywords bioregion
#' @importFrom stats quantile
#' @return a vector with the discretized values.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#' @seealso \code{\link{coldspots}}
#' @examples
#' library(sp)
#' s <- readRDS(system.file("nigeria/SR_Naija.rds", package = "phyloregion"))
#' k <- 10
#' COLOUR <- hcl.colors(k, "RdYlBu")
#' y <- choropleth(s$SR, k)
#'
#' plot(s$SR, y)
#' ## To plot and color according to some metric:
#' plot(s, col = COLOUR[y])
#' @export
#'
choropleth <- function(x, k = 10, style = "quantile", ...) {
  quants <- switch(style,
                   quantile = quantile(x, probs = seq(0, 1, 1/k), ...),
                   equal = seq(min(x), max(x), length.out = (k + 1)),
                   pretty = c(pretty(x, k = k, ...)),
                   jenks = JENKS(x, k))
  l <- length(quants) - 1
  col_vec <- integer(length(x))
  col_vec[x == quants[1]] <- 1L
  for (i in seq_len(l)) {
    col_vec[x > quants[i] & x <= quants[i + 1]] <- i
  }
  col_vec
}
