#' Map species' trait values in geographic space
#'
#' \code{map_trait} add species trait values to species distribution
#' in geographic space.
#'
#' @param x A community data object - a vector (with names matching trait data)
#' or a data.frame or matrix (with column names matching names in trait data)
#' @param trait A data.frame of species traits with a column of species names
#' matching species names in the community data, and another column with
#' the trait values.
#' @param FUN The function used to aggregate species trait values
#' in geographic space. By default, if \code{FUN = sum}, the sum of
#' all species traits per area  or grid cell is calculated.
#' @param shp a polygon shapefile of grid cells.
#' @param \dots Further arguments passed to or from other methods.
#'
#' @rdname map_trait
#' @keywords phyloregion
#' @importFrom raster values
#' @importFrom stats aggregate
#'
#' @return A data frame of species traits by site.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' data(africa)
#' x <- EDGE(africa$IUCN, africa$phylo, Redlist = "IUCN", species="Species")
#' y <- map_trait(africa$comm, x, FUN = sd, shp=africa$polys)
#'
#' plot_swatch(y, y$traits,
#'            col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))
#' @export map_trait
map_trait <- function(x, trait, FUN = sum, shp = NULL, ...){
  if(is(x, "sparseMatrix")) x <- sparse2long(x)
  grids <- NULL
  ind1 <- intersect(x$species, names(trait))
  if(length(ind1)>0){
    submat <- x[x$species %in% ind1,]
    submat <- cbind(submat, trait=trait[submat$species])
    #tmp <- data.table::as.data.table(submat)
    res <- stats::aggregate(submat$trait, by = list(submat$grids), FUN = FUN)
    #res <- aggregate(trait ~ grids, data = submat, FUN = FUN)
    #res <- tmp[, FUN(trait), by=grids]
    names(res) <- c("grids", "traits")
    if(length(shp)==0){
      m <- res
      m
    } else {
      m <- sp::merge(shp, res, by="grids")
      if (!inherits(m, "SpatialPolygons")){
        stop("Invalid geometry, may only be applied to polygons")
      }
      m <- m[!is.na(m@data$traits),]
      m
    }

  } else {
    stop("Taxa names in community data do not match names in trait data")
  }
}
