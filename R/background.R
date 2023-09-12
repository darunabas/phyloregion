#' Sample background points
#' 
#' Generates background sample as null for species distribution 
#' modeling and other things.
#' 
#' @param calib A SpatRaster of the species calibration area. If a polygon, 
#' convert to SpatRaster using rasterize.
#' @param spatkde A weighted or unweighted Gaussian Kernel Density estimate 
#' (KDE) for all input occurrence records.
#' @param size Size of background sample points
#' @rdname backg
#' @importFrom terra nlyr xyFromCell
#' @importFrom terra global crop resample mask as.data.frame
#' @return A dataframe containing the generated background points.
#' @export
backg <- function(calib, spatkde, size=10000) {
  s <- terra::global(calib, sum, na.rm=TRUE)[[1]]
  if ((s) > size) {
    kde_crop <- terra::crop(spatkde, calib)
    kde_crop <- terra::resample(kde_crop, calib)
    kde_mask <- terra::mask(kde_crop, calib)
    bg <- randpoints(kde_mask, size, prob = TRUE)
  } else {
    bg <- as.data.frame(calib, xy=TRUE)[, 1:2]
  }
  bg <- as.data.frame(bg)
  return(bg)
}

#' Random points with probability weights
#' 
#' Sample random background points using a vector of probability 
#' weights.
#' 
#' @param ras Input SpatRaster
#' @param size An positive integer of the number of samples to generate.
#' @param prob Vector of probability weights for obtaining the points sampled.
#' @rdname backg
#' @importFrom terra nlyr xyFromCell
#' @return A dataframe with sampled points.
#' @export
randpoints <- function(ras, size, prob = NULL) {
  if(terra::nlyr(ras) > 1) ras <- ras[[1]]
  v <- as.vector(ras)
  v.notNA <- which(!is.na(v))
  v.vals <- v[v.notNA]
  if(prob == TRUE & size <= length(v.notNA)) {
    res <- sample(v.notNA, size, prob = v.vals)
  } else {
    res <- sample(v.notNA, size)
  }
  terra::xyFromCell(ras, res)
}

#' Get current directory
#' 
#' Gets the path of the current directory.
#' 
#' @param path Character vector of the directory path names.
#' @rdname dirpath
#' @return A character vector containing the name of the current directory.
#'
#' @examples
#' # Get the name of the current working directory
#' dirpath()
#'
#' # Get the name of a specific directory from a path
#' dirpath("/path/to/directory")
#'
#' @export
dirpath <- function(path){
  res <- strsplit(path, "\\/")[[1]]
  res[[length(res)]]
}

