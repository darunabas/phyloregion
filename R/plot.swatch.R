#' Plot shapefile polygons based on slot values
#'
#' \code{plot_swatch} maps discretized values of a quantity using
#' continuous color gradients.
#'
#' @param x A data frame or object of the class SpatialPolygonsDataFrame
#' @param values Variable in the SpatialPolygonsDataFrame for which to
#' discretize the values of the quantity.
#' @param col A vector of colors
#' @param leg Numeric, length of the legend
#' @param key_label label for the color key
#' @param legend logical indicating whether to add a legend to the map.
#' @param breaks one of \dQuote{equal}, \dQuote{pretty}, \dQuote{jenks},
#' \dQuote{quantile} or numeric with the actual breaks by specifying
#' the minimum (\code{min}) and maximum (\code{max}) bounds.
#' @param pos location to position the legend such as \dQuote{bottomright},
#' \dQuote{bottomleft}, \dQuote{topleft}, and \dQuote{topright}.
#' @param lwd numeric, line width of the legend.
#' @param border plot polygons in SpatialPolygons object
#' @param min the minima of the lowest bound of the break.
#' @param max the maxima of the upper bound of the break
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot_swatch
#' @keywords visualization and mapping
#' @importFrom stats quantile
#' @importFrom graphics legend par points rect segments strheight strwidth text
#' @importFrom graphics xinch yinch
#' @importFrom grDevices rgb hcl.colors as.graphicsAnnot xy.coords
#' @importFrom stats median
#' @return Returns no value, just map swatch of colors in geographic space!
#' @seealso \code{\link[sp]{SpatialPolygons-class}}
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' library(sp)
#' s <- readRDS(system.file("nigeria/SR_Naija.rds", package = "phyloregion"))
#' plot_swatch(s, values = s$SR)
#' @export
plot_swatch <- function(x, values=NULL,
                        col = hcl.colors(n=10, palette = "Blue-Red 3", rev=FALSE),
                        key_label = "", leg = 10,
                        lwd = 15, pos = "bottomleft", legend = TRUE,
                        border=par("fg"), breaks = "quantile",
                        min = NULL, max = NULL, ...) {
  ed <- FALSE
  k <- length(col)
  #rev <- FALSE
  if(inherits(x, "phyloregion")){
    x <- x$evol_distinct
    if(is.null(values)) values <- x$ED
    k <- nrow(x)
    ed <- TRUE
    rev <- TRUE
  }
  if(is.null(values)) stop("You need to supply value argument!")
  x$values <- values
  #colrs <- hcl.colors(k, palette = palette, rev=rev)
  y <- choropleth(values, k, breaks = breaks, min, max) # , style = style
  plot(x, col = col[y], border = border, ...)
  if(ed) text(x, labels = as.character(x@data$cluster), ...)
  if (legend) {
    color_key(x, col, vals = values, lab = key_label, pos = pos, leg = leg,
              lwd = lwd)
  }
}

