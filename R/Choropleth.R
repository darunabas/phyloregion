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
choropleth <- function(x, k = 10, style="quantile", ...) {
  quants <- switch(style,
                   quantile = quantile(x, seq(0, 1, length.out = k + 1),
                   equal = seq(min(x), max(x), length.out = k),
                   pretty = c(pretty(x, k = k + 1))))
  l <- length(quants) - 1
  col_vec <- integer(length(x))
  col_vec[x == quants[1]] <- 1L
  for (i in seq_len(l)) {
    col_vec[x > quants[i] & x <= quants[i + 1]] <- i
  }
  col_vec
}
