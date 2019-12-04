#' Computes biodiversity hotspots
#'
#' \code{hotspots} map areas or grid cells with highest values for a biodiversity metric e.g.
#' species richness, species endemism or degree of threat.
#'
#' @param x A data frame
#' @param values Variable in the dataframe on which to compute hotspots analysis
#' @param prob The threshold quantile for representing the highest proportion of biodiversity
#' in an area. By default, the threshold is set to \code{prob = 2.5} percent.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname hotspots
#' @keywords bioregion
#' @importFrom stats quantile
#'
#' @return Integers of 1s and 0s with 1 corresponding to the hotspots.
#'
#' @references
#' Myers, M., Mittermeier, R.A., Mittermeier, C.G., da Fonseca, G.A.B. & Kent, J.
#' (2000) Biodiversity hotspots for conservation priorities. \emph{Nature}
#' \strong{403}: 853–858.
#'
#' Ceballos, G. & Ehrlich, P.R. (2006) Global mammal distributions, biodiversity
#' hotspots, and conservation. \emph{Proceedings of the National Academy of
#' Sciences USA} \strong{103}: 19374–19379.
#'
#' Orme, C.D., Davies, R.G., Burgess, M., Eigenbrod, F., Pickup, N. et al. (2005)
#' Global hotspots of species richness are not congruent with endemism or
#' threat. \emph{Nature} \strong{436}: 1016–1019.
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#' @seealso \code{\link{coldspots}}
#'
#' @examples
#' require(raster)
#' data(africa)
#' names(africa)
#'
#' Endm <- weighted.endemism(africa$comm)
#'
#' H <- hotspots(Endm, values = Endm$WE)

## Merge endemism values to shapefile of grid cells.
#' m <- merge(africa$polys, H, by="grids")
#' m <- m[!is.na(m@data$values),]
#'
#' par(mfrow = c(1,2))
#' par(mar = rep(0, 4))
#' plot(africa$polys, border="grey", col="lightgrey", main="Endemism Hotspots")
#' plot(m[(m@data$hot==1),], col="red", add=TRUE, border=NA)
#'
#' plot.swatch(m, m$WE, k=20, pos = "bottomright")

#' @export hotspots
hotspots <- function(x, values, prob = 2.5, ...){
  quant <- (1-(prob/100))
  x$values <- values
  r <-quantile(values, quant, na.rm=TRUE)
  values[which(values < r[[1]])] <- 0
  values[which(values > r[[1]])] <- 1
  values[which(values == r[[1]])] <- 1
  x$hot <- values
  x
}

