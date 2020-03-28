color_key <- function(y, cols, vals, lab = "ED", leg = 5, lwd = 15,
                      pos = "bottomright") {
  if (pos == "bottomright") {
    a <- y@bbox[3] + 0.5
    b <- y@bbox[2]
  }
  if (pos == "topleft") {
    a <- y@bbox[1] - 0.5
    b <- y@bbox[4] - leg
  }
  if (pos == "bottomleft") {
    a <- y@bbox[1] - 0.5
    b <- y@bbox[2]
  }
  if (pos == "topright") {
    a <- y@bbox[3] + 0.5
    b <- y@bbox[4] - leg
  }
  l <- length(cols)
  x <- cbind(rep(a, l), rep(a, l))
  y <- b + cbind(0:(l - 1) / l, seq_len(l) / l) * (leg)
  for (i in seq_len(l)) lines(x[i, ], y[i, ], col = cols[i], lwd = lwd,
                              lend = 2)
  text(x = a, y = b, round(min(vals), 3), pos = 4, cex = 0.7) # lim texts
  text(x = a, y = b + (leg / 2), round(median(vals), 3), pos = 4, cex = 0.7)
  text(x = a, y = b + leg, round(max(vals), 3), pos = 4, cex = 0.7)
  text(x = a, y = b + leg, lab, pos = 3, cex = 1)
}

#' Visualize biogeographic patterns
#'
#' @param x an object of class phyloregion from \code{phyloregion}
#' @param palette name of the palette to generate colors from.
#' The name is matched to the list of available color palettes from
#' the \code{hcl.colors} function in the \code{grDevices} package.
#' @param key_label label for the color key
#' @param legend logical indicating whether to add a legend to the map.
#' @param pos location to position the legend such as \dQuote{bottomright},
#' \dQuote{bottomleft}, \dQuote{topleft}, and \dQuote{topright}.
#' @param leg parameter of the color key.
#' @param lwd parameter of the color key.
#' @param shp a polygon shapefile of grid cells.
#' @param \dots arguments passed among methods.
#' @return No return value, called for plotting.
#' @rdname plot_evoldistinct
#' @export
plot_evoldistinct <- function(x, palette = "YlOrBr", pos = "bottomleft",
                              key_label = "", legend = TRUE, leg = 5,
                              lwd = 15, ...) {
  if (!inherits(x, "phyloregion"))
    stop("object \"x\" is not of class \"phyloregion\"")
  m1 <- x$shp
  k <- nrow(m1)
  COLRS <- hcl.colors(k, palette, rev = TRUE, ...)
  y <- choropleth(m1$ED, k) # , style = style
  plot(m1, col = COLRS[y], ...)
  text(m1, labels = as.character(m1@data$cluster), ...)
  if (legend) {
    color_key(m1, COLRS, vals = m1$ED, leg = leg, lwd = lwd, pos = pos,
      lab = key_label)
  }
}


#' @rdname plot_evoldistinct
#' @importMethodsFrom sp plot
#' @importMethodsFrom raster text
#' @method plot phyloregion
#' @importFrom graphics plot plot.default
#' @export
plot.phyloregion <- function(x, shp=NULL, palette="YlOrBr", ...) {
  if (!inherits(x, "phyloregion"))
    stop("object \"x\" is not of class \"phyloregion\"")
  if(is.null(shp)) shp <- x$shp
  if(is.null(shp)) stop("Need shape file!")
  k <- x$k
  if(palette=="NMDS") clr <- shp$COLOURS
    else clr <- hcl.colors(k, palette)[shp$cluster]
  plot(shp, col = clr, ...)
  text(shp, labels = as.character(shp@data$cluster), ...)
}

#' @rdname plot_evoldistinct
#' @export
plot_NMDS <- function(x, ...) {
  if (!inherits(x, "phyloregion"))
    stop("object \"x\" is not of class \"phyloregion\"")
  c1 <- x$NMDS
  plot(c1$points, bg = hexcols(c1), pch = 21, ...)
}

#' @rdname plot_evoldistinct
#' @examples
#' data(africa)
#' tree <- africa$phylo
#' x <- africa$comm
#'
#' subphy <- match_phylo_comm(tree, x)$phy
#' submat <- match_phylo_comm(tree, x)$com
#'
#' pbc <- phylobeta(submat, subphy)
#' y <- phyloregion(pbc[[1]], shp=africa$polys)
#'
#' plot_NMDS(y, cex=6)
#' text_NMDS(y, cex=2)
#' plot(y, cex=1, palette="NMDS")
#' plot(y, cex=1)
#' plot_evoldistinct(y)
#' @export
text_NMDS <- function(x, ...) {
  if (!inherits(x, "phyloregion"))
    stop("object \"x\" is not of class \"phyloregion\"")
  c1 <- x$NMDS
  text(c1$points, as.character(1:c1$nobj), ...)
}
