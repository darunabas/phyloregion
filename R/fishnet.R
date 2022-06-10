#' Create a fishnet of regular grids
#'
#' The \code{fishnet} function creates a regular grid of locations covering
#' the study area at various grain sizes.
#'
#' @param mask a polygon shapefile covering the boundary of the survey region.
#' @param res the grain size of the grid cells in decimal degrees (default).
#' @param type the type of grid cell; either \dQuote{square} (default) for
#' square grids, or \dQuote{hexagon} for hexagonal grids.
#' @rdname fishnet
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values extent res<-
#' @importFrom sp CRS proj4string spsample HexPoints2SpatialPolygons
#' @importFrom sp SpatialPolygonsDataFrame
#' @return A spatial polygon object of equal area grid cells covering the
#' defined area.
#' @references
#' Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum entropy
#' modeling of species geographic distributions. \emph{Ecological Modelling}
#' \strong{190}: 231-259.
#' @examples
#' file <- system.file("nigeria/nigeria.rds", package="phyloregion")
#' d <- readRDS(file)
#' d1 <- fishnet(d, res = 0.75)
#' @export
fishnet <- function(mask, res = 0.5, type = "square"){

  res <- switch(type,
                square = {
                  s <- raster(extent(mask))
                  res(s) <- res
                  suppressWarnings(invisible(proj4string(s) <- proj4string(mask)))
                  m <- rasterToPolygons(s)
                  m$grids <- paste0("v", seq_len(nrow(m)))
                  m <- m[, "grids"]
                  spo <- m[subset(mask), ]
                },
                hexagon = {
                  r <- suppressWarnings(invisible(spsample(mask,
                                                           type = "hexagonal",
                                                           cellsize = res)))
                  f <- HexPoints2SpatialPolygons(r)

                  spo <- SpatialPolygonsDataFrame(f,
                        data.frame(grids=paste0("v", seq_along(f))),
                        match.ID = FALSE)
                })
  return(res)
}
