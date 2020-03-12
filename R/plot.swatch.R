#' Plot shapefile polygons based on slot values
#'
#' \code{plot_swatch} maps discretized values of a quantity based on their
#' quantiles.
#'
#' @param x A data frame or object of the class SpatialPolygonsDataFrame
#' @param values Variable in the SpatialPolygonsDataFrame for which to
#' discretize the values of the quantity.
#' @param k Numeric, the desired number of bins to discretize.
#' @param palette name of the palette to generate colors from. The name is
#' matched to the list of available color palettes from the \code{hcl.colors}
#' function in the \code{grDevices} package.
#' @param leg Numeric, length of the legend
#' @param key_label label for the color key
#' @param legend logical indicating whether to add a legend to the map.
#' @param pos location to position the legend such as \dQuote{bottomright},
#' \dQuote{bottomleft}, \dQuote{topleft}, and \dQuote{topright}.
#'
#' @param lwd numeric, line width of the legend.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot_swatch
#' @keywords visualization and mapping
#' @importFrom stats quantile
#' @importFrom graphics legend par points rect segments strheight strwidth text
#' @importFrom graphics xinch yinch
#' @importFrom sp plot
#' @importFrom grDevices rgb hcl.colors as.graphicsAnnot xy.coords
#' @importFrom stats median
#' @return Returns no value, just map swatch of colors in geographic space!
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' library(raster)
#' s <- readRDS(system.file("nigeria/SR_Naija.rds", package = "phyloregion"))
#' plot_swatch(s, values = s$SR, k = 20)
#' @export
plot_swatch <- function(x, values, k = 10, palette = "Blue-Red 3",
                        key_label = "", leg = 10, lwd = 15, pos = "bottomleft",
                        legend = TRUE, ...) {
  x$values <- values
  colrs <- hcl.colors(k, palette = palette, ...)
  y <- choropleth(values, k, ...) # , style = style
  raster::plot(x, col = colrs[y], border = NA, ...)
  if (legend) {
    color_key(x, colrs, vals = values, lab = key_label, pos = pos, leg = leg,
              lwd = lwd)
  }
}
