#' Overlay of polygons to community matrix
#'
#' This function performs spatial overlays of polygons with polygons
#' to generate a community matrix
#'
#' @param dat layers of merged maps corresponding to species ranges from
#' which the geometries or attributes are queried.
#' @param mask a polygon shapefile covering the boundary of the survey region.
#' @param res the grain size of the grid cells in degrees.
#' @param shp.grids if specified, the polygon shapefile of grid cells
#' @param species a character string. The column with the species or taxon name.
#' Default = \dQuote{species}.
#' with a column labeled "grids"
#' @rdname polys2comm
#'
#' @param \dots Further arguments passed to or from other methods.
#' @keywords phyloregion
#' @importFrom data.table rbindlist
#' @importFrom sp coordinates over CRS proj4string merge
#'
#' @inheritParams fishnet
#'
#' @return
#' \itemize{
#'   \item comm_dat: community data frame
#'   \item poly_shp: shapefile of grid cells with the values per cell.
#' }
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' require(raster)
# some example data
#' p1 <- Polygons(list(Polygon(rbind(c(-180,-20), c(-140,55), c(10, 0),
#'                                   c(-140,-60), c(-180,-20)))), 1)
#' p2 <- Polygons(list(Polygon(rbind(c(-10,0), c(140,60), c(160,0),
#'                                   c(140,-55), c(-10,0)))), 2)
#' p3 <- Polygons(list(Polygon(rbind(c(-125,0), c(0,60), c(40,5),
#'                                   c(15,-45), c(-125,0)))), 3)
#' sp <- SpatialPolygons( list( p1 , p2, p3),
#'                        proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
#' sp <- SpatialPolygonsDataFrame(sp, data.frame(Species=c("sp1", "sp2", "sp3")))
#'
#' pol <- polys2comm(dat = sp, species="Species")
#' plot_swatch(pol$poly_shp, values = pol$poly_shp$richness, k=10)
#' @export
polys2comm <- function(dat,
                       res=1,
                       shp.grids=NULL,
                       species = "species",
                       verbose = TRUE, ...){
  if (verbose) {
    message("Generating community data from polygon shapefiles")
  }

  dat <- dat[, species, drop=FALSE]
  names(dat) <- "species"
  if(length(shp.grids)==0){

    e <- extent(dat)+(2*res)
    # coerce to a SpatialPolygons object
    mask <- as(e, 'SpatialPolygons')
    lu <- as.data.frame(1L)
    mask <- sp::SpatialPolygonsDataFrame(mask, lu)
    m <- fishnet(mask, res = res)
  } else (m = shp.grids)
  proj4string(dat) <- proj4string(m)
  x <- mapply(cbind, sp::over(dat, m, returnList = TRUE),
              dat@data$species, SIMPLIFY=F)
  y <- rbindlist(x)
  names(y) <- c("grids", "species")
  y <- y[complete.cases(y),]
  y <- unique(y[, c("grids", "species")])
  res <- data.frame(table(y$grids))
  names(res) <- c("grids", "richness")
  z <- sp::merge(m, res, by="grids")
  z <- z[!is.na(z@data$richness),]
  return(list(comm_dat=y, poly_shp=z, grids=m))
}

