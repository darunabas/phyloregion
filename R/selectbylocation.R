#' Select polygon features from another layer and adds polygon attributes to layer
#'
#' The \code{selectbylocation} function selects features based on
#' their location relative to features in another layer.
#'
#' @param x source layer of the class SpatVect
#' @param y Target layer or mask extent to subset from.
#' @rdname selectbylocation
#' @keywords bioregion
#' @importFrom terra crs<- project crs
#'
#' @export
#' @return A spatial polygons or spatial points object pruned to the extent
#' of the target layer.
#'
#' @examples
#' library(terra)
#' d <- vect(system.file("ex/nigeria.json", package="phyloregion"))
#' e <- ext(d)
#'
#' set.seed(1)
#' m <- data.frame(lon = runif(1000, e[1], e[2]),
#'                 lat = runif(1000, e[3], e[4]),
#'                 sites = seq(1000))
#' m <- vect(m)
#' z <- selectbylocation(m, d)
#' plot(d)
#' points(m, col = "blue", pch = "+")
#' points(z, col = "red", pch = "+")
selectbylocation <- function(x, y) {
  pj <- "+proj=longlat +datum=WGS84"
  crs(x) <- pj
  crs(y) <- crs(x)
  y <- project(y, pj)
  x[y, ]
}
