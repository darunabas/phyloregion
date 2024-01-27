#' Standardizes raster values for mapping
#'
#' This function standardizes values of a raster layer for mapping.
#'
#' @param ras an input raster layer.
#' @param ref a raster layer for reference.
#' @importFrom terra as.data.frame setValues
#' @return A raster layer that has been standardized and ready for mapping 
#
#' @export
rast_quantile <- function(ras, ref) {
  x <- as.data.frame(ras, na.rm=FALSE)
  x$grids <- paste0("v", seq_len(nrow(x)))
  x <- na.omit(x)
  nm <- names(x)[1]
  y <- (rank(x[[nm]])/length(x[[nm]]))*100
  names(y) <- x$grids
  q <- ras
  values(q) <- paste0("v", seq_len(ncell(q)))
  i <- match(as.data.frame(q)[[1]], names(y))
  z <- setValues(q, y[i]) |> crop(ref) |> mask(ref)
  return(z)
}

