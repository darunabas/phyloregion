#' Bin shapefile polygons based on slot values
#'
#' \code{choropleth} discretizes the values of a quantity based on their quantiles. 
#'
#' @param x A data frame or object of the class SpatialPolygonsDataFrame
#' @param values Variable in the SpatialPolygonsDataFrame for which to discretize 
#' the values of the quantity.
#' @param k Numeric, the desired number of bins to discretize.
#' @param style one of \dQuote{equal}, \dQuote{pretty}, or \dQuote{quantile}.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname choropleth
#' @keywords bioregion
#' @importFrom stats quantile 
#'
#' @export
#' @return
#' {returns a SpatialPolygonsDataFrame with a column of the discretized values}
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' require(raster)
#' s <- readRDS(system.file("nigeria/SR_Naija.rds", package= "phyloregion"))
#' k=10
#' COLOUR <- hcl.colors(k, "RdYlBu") 
#' y = choropleth(s, values=s$SR, k)
#' 
#' ## To plot and color according to some metric:
#' plot(y, col=COLOUR[y$values], border = NA)
#'
#' 
choropleth <- function(x, values, k=10, ...){
  x$values <- values
  quants <- quantile(values, seq(0,1, length.out = k+1))
  l = length(quants) -1
  col_vec <- numeric(length(values))
  col_vec[values==quants[1]] <- 1
  for(i in seq_len(l)){
    col_vec[values>quants[i] & values<= quants[i+1] ] = i
  }
  x$values <- col_vec
  x
}
