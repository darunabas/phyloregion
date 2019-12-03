legend <- function(y,
                   cols,
                   vals,
                   lab="ED",
                   leg=5,
                   lwd=15,
                   pos="bottomright") {

  if(pos=="bottomright"){
    a=y@bbox[3]+0.5
    b=y@bbox[2]
  }
  if(pos=="topleft"){
    a=y@bbox[1]-0.5
    b=y@bbox[4]-leg
  }
  if(pos=="bottomleft"){
    a=y@bbox[1]-0.5
    b=y@bbox[2]
  }
  if(pos=="topright"){
    a=y@bbox[3]+0.5
    b=y@bbox[4]-leg
  }
  X <- cbind(rep(a, length(cols)), rep(a, length(cols)))
  Y <- b + cbind(0:(length(cols) - 1)/length(cols), 1:length(cols)/length(cols)) *
    (leg)
  for (i in 1:length(cols)) lines(X[i, ], Y[i, ], col = cols[i],
                                  lwd = lwd, lend = 2)
  text(x = a, y = b, round(min(vals), 3), pos = 4,
       cex = 0.7) # lim texts
  text(x = a, y = b+(leg/2), round(median(vals), 3), pos = 4,
       cex = 0.7)
  text(x = a, y = b+leg, round(max(vals), 3),
       pos = 4, cex = 0.7)
  text(x = a, y = b+leg, lab,
       pos = 3, cex = 1)
}


#' Plot shapefile polygons based on slot values
#'
#' \code{plot_swatch} maps discretized values of a quantity based on their quantiles.
#'
#' @param x A data frame or object of the class SpatialPolygonsDataFrame
#' @param values Variable in the SpatialPolygonsDataFrame for which to discretize
#' the values of the quantity.
#' @param k Numeric, the desired number of bins to discretize.
#' @param swatch name of the palette to generate colors from. The name is matched
#' to the list of available color palettes from the code{hcl.colors} function in
#' the code{grDevices} package.
#' @param leg Numeric, length of the legend
#' @param legend logical indicating whether to add a legend to the map.
#' @param pos location to position the legend such as \dQuote{bottomright},
#' \dQuote{bottomleft}, \dQuote{topleft}, and \dQuote{topright}.
#'
#' @param lwd numeric, line width of the legend.
#' @param title a character string indicating the caption to be placed on top
#' of the legend.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot_swatch
#' @keywords visualization and mapping
#' @importFrom stats quantile
#' @importFrom graphics legend par points rect segments strheight strwidth text xinch yinch
#' @importFrom sp plot
#' @importFrom grDevices rgb hcl.colors as.graphicsAnnot xy.coords
#' @importFrom stats median
#'
#' @inheritParams choropleth
#'
#' @return Returns no value, just map swatch of colors in geographic space!
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' require(raster)
#' s <- readRDS(system.file("nigeria/SR_Naija.rds", package= "phyloregion"))
#' plot_swatch(s, values = s$SR, k=20)
#' @export
plot_swatch <- function(x, values, k = 10, swatch = "Blue-Red 3", lab="",
                         leg=5, lwd=15, pos="bottomright", legend=TRUE, ...)
{
  x$values <- values
  COLRS <- hcl.colors(k, swatch)
  y <- choropleth(x, values, k) #, style = style
  plot(y, col = COLRS[y$values], border = NA)
  if(legend) {
    legend(y, COLRS, vals=values, lab=lab, pos=pos, leg=leg, lwd=lwd)
  }
}
