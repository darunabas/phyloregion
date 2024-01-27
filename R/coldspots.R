#' Computes biodiversity coldspots and hotspots
#'
#' \code{coldspots} and \code{hotspots} map areas or grid cells with lowest
#' or highest values, respectively, of a biodiversity metric e.g.
#' species richness, species endemism or degree of threat.
#'
#' @param x a vector on which to compute hotspots or coldspots
#' @param y a vector on which to compare x against
#' @param prob The threshold quantile for representing the lowest
#' (\code{coldspots}) or highest (\code{hotspots}) proportion of biodiversity in
#' an area. By default, the threshold is set to \code{prob = 2.5} percent.
#' @param ras a SpatRaster on which to compute hotspots.
#' @param ref a raster layer for reference.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname coldspots
#' @keywords phyloregion
#' @importFrom stats quantile
#' @importFrom terra as.data.frame setValues
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
#'
#' @examples
#' library(terra)
#' data(africa)
#' p <- vect(system.file("ex/sa.json", package = "phyloregion"))
#'
#' Endm <- weighted_endemism(africa$comm)
#' C <- coldspots(Endm, na.rm=TRUE) # coldspots
#' H <- hotspots(Endm, na.rm=TRUE) # hotspots
#'
#' ## Merge endemism values to shapefile of grid cells.
#' DF <- data.frame(grids=names(C), cold=C, hot=H)
#' m <- merge(p, DF, by = "grids", all = TRUE)
#'
#' plot(p, border = "grey", col = "lightgrey",
#'      main = "Weighted Endemism Hotspots and Coldspots")
#' plot(m[(m$cold == 1), ], col = "blue", add = TRUE, border = NA)
#' plot(m[(m$hot == 1), ], col = "red", add = TRUE, border = NA)
#' legend("bottomleft", fill = c("blue", "red", "yellow", "green"),
#'        legend = c("coldspots", "hotspots"), bty = "n", inset = .092)
coldspots <- function(x, y = NULL, prob = 2.5, ...) {
  quant <- prob / 100
  if(!is.null(y)) {
    r <- quantile(y, quant, ...)
  } else {
    r <- quantile(x, quant, ...)
  }
  x[which(x < r[[1]])] <- NA
  x[which(x > r[[1]])] <- 0
  x[which(x == r[[1]])] <- NA
  x[which(is.na(x))] <- 1
  x
}


#' @rdname coldspots
#' @export
hotspots <- function(x, y = NULL, prob = 2.5, ...) {
  quant <- (1 - (prob/100))
  if(!is.null(y)) {
    r <- quantile(y, quant, ...)
  } else {
    r <- quantile(x, quant, ...)
  }
  x[which(x < r[[1]])] <- NA
  x[which(x > r[[1]])] <- 1
  x[which(x == r[[1]])] <- 1
  x[which(is.na(x))] <- 0
  x
}

#' @rdname coldspots
#' @export
rast_hotspot <- function(ras, ref=NULL, prob=10) {
  x <- as.data.frame(ras, na.rm=FALSE)
  x$grids <- paste0("v", seq_len(nrow(x)))
  nm <- names(ras)
  y <- hotspots(x=x[[nm]], prob = prob, na.rm=TRUE)
  names(y) <- x$grids
  q <- ras
  values(q) <- paste0("v", seq_len(ncell(q)))
  i <- match(as.data.frame(q)[[1]], names(y))
  if(!is.null(ref)) {
    z <- setValues(q, y[i]) |> crop(ref) |> mask(ref)
  } else {
    z <- setValues(q, y[i]) |> crop(ras) |> mask(ras)
  }
  names(z) <- nm
  return(z)
}

