#' Computes biodiversity coldspots
#'
#' \code{coldspots} map areas or grid cells with lowest values of a biodiversity metric e.g. 
#' species richness, species endemism or degree of threat.
#'
#' @param x A data frame
#' @param values Variable in the dataframe on which to compute coldspots analysis
#' @param prob The threshold quantile for representing the lowest proportion of biodiversity 
#' in an area. By default, the threshold is set to \code{prob = 2.5} percent.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname coldspots
#' @keywords phyloregion
#' @importFrom stats quantile
#'
#' @export
#' @return
#' {Integers of 1s and 0s with 1 corresponding to the coldspots}
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
#' @seealso \code{\link[bioregion]{hotspots}}
#'
#' @examples
#' require(raster)
#' data(africa)
#' names(africa)
#' 
#' Endm <- weighted.endemism(africa$comm)
#' C <- coldspots(Endm, values = Endm$WE) # coldspots
#' H <- hotspots(Endm, values = Endm$WE) # hotspots

## Merge endemism values to shapefile of grid cells.
#' m <- Reduce(function(x,y) sp::merge(x,y,by="grids",all=TRUE) , 
#'                           list(africa$polys, C, H))
#' m <- m[!is.na(m@data$values.x),]
#' 
#' plot(africa$polys, border="grey", col="lightgrey", main="Weighted Endemism Hotspots and Coldspots")
#' plot(m[(m@data$cold==1),], col="blue", add=TRUE, border=NA)
#' plot(m[(m@data$hot==1),], col="red", add=TRUE, border=NA)
#' legend("bottomleft", fill = c("blue", "red", "yellow", "green"), 
#'        legend = c("coldspots", "hotspots"), bty = "n", inset=.092)
#' 
coldspots <- function(x, values, prob = 2.5, ...){
  quant <- prob/100
  x$values <- values
  r <- quantile(values, quant, na.rm=TRUE)
  values[which(values < r[[1]])] <- NA
  values[which(values > r[[1]])] <- 0 
  values[which(values == r[[1]])] <- NA
  values[which(is.na(values))] <- 1
  x$cold <- values
  x
}

