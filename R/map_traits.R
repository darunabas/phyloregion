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
#' @param pol a vector polygon of grid cells.
#' @param \dots Further arguments passed to or from other methods.
#'
#' @rdname map_trait
#' @keywords phyloregion
#' @importFrom stats aggregate
#' @importFrom terra merge geomtype
#'
#' @return A data frame of species traits by site.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' data(africa)
#' library(terra)
#' p <- vect(system.file("ex/sa.json", package = "phyloregion"))
#' x <- EDGE(africa$IUCN, africa$phylo, Redlist = "IUCN",
#'           species = "Species")
#' y <- map_trait(africa$comm, x, FUN = sd, pol = p)
#'
#' plot(y, "traits", col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))
#' @export map_trait
map_trait <- function(x, trait, FUN = sum, pol = NULL, ...){
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
    if(length(pol)==0){
      m <- res
      m
    } else {
      m <- terra::merge(pol, res, by="grids")
      if (!geomtype(m)=="polygons") {
        stop("Invalid geometry, may only be applied to vector polygons")
      }
      m
    }

  } else {
    stop("Taxa names in community data do not match names in trait data")
  }
}
