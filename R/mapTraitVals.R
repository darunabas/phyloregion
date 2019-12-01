#' Map species trait values in geographic space
#'
#' \code{mapTraits} add species trait values to species distribution in geographic space.
#'
#' @param x A community data object - a vector (with names matching trait data)
#' or a data.frame or matrix (with column names matching names in trait data)
#' @param trait A data.frame of species traits with a column of species names
#' matching species names in the community data, and another column with the trait values.
#' @param FUN The function used to aggregate species trait values in geographic space.
#' By default, if \code{FUN = sum}, the sum of all species traits per area or grid cell
#' is calculated.
#' @param \dots Further arguments passed to or from other methods.
#'
#' @rdname mapTraits
#' @keywords phyloregion
#' @importFrom raster values
#' @importFrom data.table as.data.table
#'
#' @return A data frame of species traits by site.
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' require(data.table)
#' fdir <- system.file("NGAplants", package="phyloregion")
#' files <- file.path(fdir, dir(fdir))
#' dat <- data.frame(raster2comm(files))
#' dd <- data.frame(table(dat$species))
#' names(dd) <- c("species", "range")
#' trait_range <- dd$range
#' names(trait_range) <- dd$species
#' gg <- mapTraits(dat, trait = trait_range, FUN = sd)
#' @export mapTraits
mapTraits <- function(x, trait, FUN = sum, ...){
  ind1 <- intersect(x$species, names(trait))
  if(length(ind1)>0){
    submat <- x[x$species %in% ind1,]
    submat <- cbind(submat, trait=trait[submat$species])
    tmp <- data.table::as.data.table(submat)
    res <- tmp[, FUN(trait), by=grids]
  } else {
    stop("Taxa names in community data do not match names in trait data")
  }
  res
}
