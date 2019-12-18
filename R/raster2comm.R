make_poly <- function(file){
  dd <- raster(file)
  pol <- rasterToPolygons(dd, fun=NULL, dissolve=FALSE, na.rm=FALSE)
  proj4string(pol) = CRS("+proj=longlat +datum=WGS84")
  pol$grids <- paste0("v", 1:nrow(pol))
  xx <- as.data.frame(xyFromCell(dd, cell=1:ncell(dd))) #Make dataframe of all xy coordinates )
  xx$grids <- paste0("v", 1:nrow(xx))

  m <- merge(pol, xx, by="grids")
  names(m)[2] <- "SR"
  names(m)[3] <- "lon"
  names(m)[4] <- "lat"
  m
}



#' Convert raster to community matrix
#'
#' \code{raster2comm} converts a set of raster files from a species distribution
#' modeling into a community matrix.
#'
#' @param files  list of raster layer objects with the same spatial
#' extent and resolution.
#' @rdname raster2comm
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values
#' @importFrom sp CRS proj4string<-
#' @examples
#' \dontrun{
#' fdir <- system.file("NGAplants", package="phyloregion")
#' files <- file.path(fdir, dir(fdir))
#' r <- raster2comm(files)
#' }
#' @export
raster2comm <- function(files) {
  poly <- make_poly(files[1])
  tmp <- raster(files[1])
  fg <- as.data.frame(xyFromCell(tmp, cell=1:ncell(tmp)))
  ind1 <- paste(as.character(fg$x),as.character(fg$y), sep="_")
  ind2 <- paste(as.character(poly$lon),as.character(poly$lat), sep="_")
  index <- match(ind1, ind2)
  res <- NULL
  #  nam <- names(obj)
  cells <- as.character(poly$grids)
  for(i in seq_along(files)) {
    obj <- raster(files[i])
    tmp <- values(obj)
#    tmp <- which(tmp>0)
    tmp <- cells[index[tmp > 0]]
    tmp <- tmp[!is.na(tmp)]
    if(length(tmp)>0)res <- rbind(res, cbind(tmp, names(obj)))
  }
  if(length(tmp)>0)colnames(res) <- c("grids", "species")
  res
}



