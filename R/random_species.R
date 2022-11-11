#' Generate random species distributions in space
#'
#' This function generates random species distributions in geographic space as
#' extent of occurrence range polygons based on convex hulls of random points.
#'
#' @param n vector of one or more elements to choose from, or a positive
#' integer.
#' @param species the desired number of species.
#' @param pol the vector polygon of the study area for determining the
#' species distributions
#' @rdname random_species
#'
#' @param \dots Further arguments passed to or from other methods.
#' @keywords phyloregion
#' @importFrom terra spatSample vect as.data.frame crop
#' @importFrom grDevices chull
#'
#' @return A vector polygon of species' extent of occurrence ranges.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @export

random_species <- function(n, species, pol, ...) {
    x <- sample(n, species, replace = TRUE, ...)
    x[x < 4] <- 4L
    y <- vector("list", length(x))
    for (i in seq_along(x)) {
        sp <- sample(c("random", "regular"), 1)
        v <- as.data.frame(terra::spatSample(pol, x[i], method = sp, ...),
                           geom = "XY")
        v <- v[, c("x", "y")]
        v <- v[chull(v), ]
        v <- as.matrix(cbind(id=1, part=1, v[, c("x", "y")]))
        s <- vect(v, type="polygons")
        e <- terra::crop(s, pol)
        if(geomtype(e)=="polygons") {
            e$species <- paste0("species", i)
            y[[i]] <- e
        }
    }
    vect(y)
}
