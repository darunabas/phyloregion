#' Select features within polygon from another layer
#'
#' The \code{selectbylocation} function selects features based on 
#' their location relative to features in another layer.
#' 
#' @param x source layer of the class SpatialPolygonsDataFrame or SpatialPointsDataFrame
#' @param y Target layer or mask extent to subset from.
#' @rdname selectbylocation
#' @keywords bioregion
#' @importFrom sp CRS proj4string coordinates<-
#'
#' @export
#' @return A spatial polygons or spatial points object pruned to the extent 
#' of the target layer.
#'
#' @examples
#' require(raster)
#' require(maptools)
#' data("wrld_simpl")
#' set.seed(1)
#' m <- data.frame(lon = runif(2000, -180, 180), 
#'                 lat = runif(2000, -90, 90), sites=seq(2000))
#' coordinates(m) <- ~lon+lat

#' z <- selectbylocation(m, wrld_simpl)

#' plot(wrld_simpl)
#' points(m, col="blue", pch="+")
#' points(z, col="red", pch="+")

selectbylocation <- function(x, y){
  proj4string(x) = proj4string(y)
  if(class(x)=="SpatialPolygonsDataFrame"){
    p <- x[subset(y), ]
  }
  else if(class(x)=="SpatialPointsDataFrame"){
    p <- x[y,]
  }
  p
}

