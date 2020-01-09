#' Create a fishnet of regular grids
#'
#' The \code{fishnet} function creates a regular grid of locations covering
#' the study area at various grain sizes.
#'
#' @param mask a polygon shapefile covering the boundary of the survey region.
#' @param res the grain size of the grid cells in decimal degrees (default).
#' @rdname fishnet
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values extent res<-
#' @importFrom sp CRS proj4string
#' @return A spatial polygon object of equal area grid cells covering the
#' defined area.
#' @references
#' Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum entropy
#' modeling of species geographic distributions. \emph{Ecological Modelling}
#' \strong{190}: 231-259.
#' @examples
#' library(raster)
#' file <- system.file("nigeria/nigeria.rds", package="phyloregion")
#' d <- readRDS(file)
#' d1 <- fishnet(d, res = 0.75)
#' @export
fishnet <- function(mask, res = 0.5){
  s <- raster(extent(mask))
  res(s) <- res
  proj4string(s) <- proj4string(mask)
  m <- rasterToPolygons(s)
  m$grids <- paste0("v", seq_len(nrow(m)))
  m <- m[, "grids"]
  dd <- m[subset(mask), ]
  dd
}
