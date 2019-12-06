#' Convert point occurrences to community data
#'
#' \code{points2comm} converts point occurrence records (e.g. from GBIF,
#' iDigBio) to a community data sample
#'
#' @param dat A point occurrence data frame with at least three columns:
#' \itemize{
#'   \item Column 1: \code{species} (listing the taxon names)
#'   \item Column 2: \code{decimallongitude} (corresponding to decimal longitude)
#'   \item Column 3: \code{decimallatitude} (corresponding to decimal latitude)
#' }
#' @param mask a polygon shapefile covering the boundary of the survey region.
#' @param res the grain size of the grid cells in degrees.
#' @param lon character with the column name of the longitude.
#' @param lat character with the column name of the lattitude.
#' @param shp.grids if specified, the polygon shapefile of grid cells
#' with a column labeled \dQuote{grids}.
#' @param species a character string. The column with the species or taxon name.
#' Default = \dQuote{species}.
#' @rdname points2comm
#'
#' @param index Diversity index, one of \dQuote{taxon_richness} (default) or
#' \dQuote{abundance}.
#' @param \dots Further arguments passed to or from other methods.
#' @keywords phyloregion
#' @importFrom sp coordinates<- over CRS proj4string merge
#' @importFrom sp coordinates<- CRS proj4string<-
#' @importFrom stats complete.cases
## @inheritParams fishnet
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
#' s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
#'
#' set.seed(1)
#' m <- data.frame(sp::spsample(s, 10000, type="nonaligned"))
#' names(m) <- c("lon", "lat")
#' species <- paste0("sp", sample(1:1000))
#' m$taxon <- sample(species, size = nrow(m), replace = TRUE)
#' pts <- points2comm(dat = m, mask = s, res = 0.5, lon="lon", lat="lat", species="taxon")
#'
#' plot_swatch(pts$poly_shp, values = pts$poly_shp$richness, k=10)
#' @export
points2comm <- function(dat, mask, res=1, lon = "decimallongitude",
                        lat = "decimallatitude", species = "species",
                        shp.grids=NULL, index="taxon_richness", ...){
  dat <- as.data.frame(dat)
  dat <- dat[, c(species, lon, lat)]
  names(dat) <- c("species", "lon", "lat")

  dat <- dat[complete.cases(dat),]

  coordinates(dat) = ~lon+lat

  if(length(shp.grids)==0){
    m <- fishnet(mask = mask, res = res)
  } else (m = shp.grids)
  proj4string(dat) <- proj4string(m)
  x <- over(dat, m)
  y <- cbind(as.data.frame(dat), x)
  y <- y[complete.cases(y),]
  if (index=="abundance"){
    y <- y[, c("grids", "species")]
  }
  else if (index == "taxon_richness") {
    y <- unique(y <- y[, c("grids", "species")])
  }
  res <- data.frame(table(y$grids))
  names(res) <- c("grids", "richness")
  z <- sp::merge(m, res, by="grids")
  z <- z[!is.na(z@data$richness),]
  return(list(comm_dat=y, poly_shp=z))
}

