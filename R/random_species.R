#' Generate random species distributions in space
#'
#' This function generates random species distributions in geographic space
#' as extent of occurrence range polygons based on convex hulls of random points.
#'
#' @param n vector of one or more elements to choose from, or a positive integer.
#' @param species the desired number of species.
#' @param shp the polygon shapefile of the study area for determining the
#' species distributions
#' @rdname random_species
#'
#' @param \dots Further arguments passed to or from other methods.
#' @keywords phyloregion
#' @importFrom rgeos gConvexHull
#' @importFrom sp spsample SpatialPolygonsDataFrame
#' @importFrom raster bind
#'
#' @return A polygon shapefile of species' extent of occurence ranges.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @export

random_species <- function(n, species, shp, ...){
    x <- sample(1:n, species, replace=TRUE, ...)

    y <- list()
    for (i in 1:length(x)) {
        sp <- sample(c("random", "regular", "stratified",
                       "nonaligned", "hexagonal", "clustered"), 6)
        v <- sp::spsample(shp, x[i], type = sp[1], iter = 10, ...)
        k <- rgeos::gConvexHull(v)
        y[[i]] <- k
    }

    m <- do.call(raster::bind, y)

    sp <- sp::SpatialPolygonsDataFrame(m,
                                   data.frame(species=paste0("species", 1:length(x))))
    sp
}

