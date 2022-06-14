#' Select polygon features from another layer and adds polygon attributes to layer
#'
#' The \code{selectbylocation} function selects features based on
#' their location relative to features in another layer.
#'
#' @param x source layer of the class SpatialPolygonsDataFrame or
#' SpatialPointsDataFrame
#' @param y Target layer or mask extent to subset from.
#' @rdname selectbylocation
#' @keywords bioregion
#' @importFrom sp CRS proj4string coordinates<-
#'
#' @export
#' @return A spatial polygons or spatial points object pruned to the extent
#' of the target layer.
#'
#' @examples
#' library(raster)
#' file <- system.file("nigeria/nigeria.rds", package = "phyloregion")
#' d <- readRDS(file)
#' e <- extent(d)
#'
#' set.seed(1)
#' m <- data.frame(lon = runif(1000, e[1], e[2]),
#'   lat = runif(1000, e[3], e[4]),
#'   sites = seq(1000))
#' coordinates(m) <- ~ lon + lat
#' z <- selectbylocation(m, d)
#' plot(d)
#' points(m, col = "blue", pch = "+")
#' points(z, col = "red", pch = "+")
selectbylocation <- function(x, y) {
  proj4string(x) <- proj4string(y)
  if (inherits(x, "SpatialPolygonsDataFrame")) {
    p <- x[subset(y), ]
  }
  else if (inherits(x, "SpatialPointsDataFrame")) {
    p <- x[y, ]
  }
  p
}
