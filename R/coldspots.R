#' Computes biodiversity coldspots and hotspots
#'
#' \code{coldspots} and \code{hotspots} map areas or grid cells with lowest
#' or highest values, respectively, of a biodiversity metric e.g.
#' species richness, species endemism or degree of threat.
#'
#' @param x a vector on which to compute coldspots
#' @param prob The threshold quantile for representing the lowest
#' (\code{coldspots}) or highest (\code{hotspots}) proportion of biodiversity in
#' an area. By default, the threshold is set to \code{prob = 2.5} percent.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname coldspots
#' @keywords phyloregion
#' @importFrom stats quantile
#' @export
#' @return A vector of integers of 1s and 0s with 1 corresponding to the
#' coldspots or hotspots
#' @references
#' Myers, M., Mittermeier, R.A., Mittermeier, C.G., da Fonseca, G.A.B. &
#' Kent, J. (2000) Biodiversity hotspots for conservation priorities.
#' \emph{Nature} \strong{403}: 853–858.
#'
#' Ceballos, G. & Ehrlich, P.R. (2006) Global mammal distributions, biodiversity
#' hotspots, and conservation. \emph{Proceedings of the National Academy of
#' Sciences USA} \strong{103}: 19374–19379.
#'
#' Orme, C.D., Davies, R.G., Burgess, M., Eigenbrod, F., Pickup, N. et al.
#' (2005) Global hotspots of species richness are not congruent with endemism or
#' threat. \emph{Nature} \strong{436}: 1016–1019.
#'
#' Daru, B.H., Van der Bank, M. & Davies, T.J. (2015) Spatial incongruence among
#' hotspots and complementary areas of tree diversity in southern Africa.
#' \emph{Diversity and Distributions} \strong{21}: 769-780.
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#' @seealso \code{\link{choropleth}}
#'
#' @examples
#' library(raster)
#' library(sp)
#' data(africa)
#' names(africa)
#'
#' Endm <- weighted_endemism(africa$comm)
#' C <- coldspots(Endm) # coldspots
#' H <- hotspots(Endm) # hotspots
#'
#' ## Merge endemism values to shapefile of grid cells.
#' DF <- data.frame(grids=names(C), cold=C, hot=H)
#' m <- merge(africa$polys, DF, by = "grids", all = TRUE)
## m <- m[!is.na(m@data$values.x), ]
#'
#' plot(africa$polys, border = "grey", col = "lightgrey",
#'   main = "Weighted Endemism Hotspots and Coldspots")
#' plot(m[(m@data$cold == 1), ], col = "blue", add = TRUE, border = NA)
#' plot(m[(m@data$hot == 1), ], col = "red", add = TRUE, border = NA)
#' legend("bottomleft", fill = c("blue", "red", "yellow", "green"),
#'   legend = c("coldspots", "hotspots"), bty = "n", inset = .092)
coldspots <- function(x, prob = 2.5, ...) {
  quant <- prob / 100
  r <- quantile(x, quant, na.rm = TRUE)
  x[which(x < r[[1]])] <- NA
  x[which(x > r[[1]])] <- 0
  x[which(x == r[[1]])] <- NA
  x[which(is.na(x))] <- 1
  x
}


#' @rdname coldspots
#' @export
hotspots <- function(x, prob = 2.5,...) {
  quant <- (1 - (prob / 100))
  r <- quantile(x, quant, na.rm = TRUE)
  x[which(x < r[[1]])] <- NA
  x[which(x > r[[1]])] <- 1
  x[which(x == r[[1]])] <- 1
  x[which(is.na(x))] <- 0
  x
}
