#' Create fishnet of regular grids
#'
#' The \code{fishnet} function creates a regular grid of locations covering
#' the study area at various grain sizes.
#'
#' @param mask a vector polygon covering the boundary of the survey region.
#' @param res the grain size of the grid cells in decimal degrees (default).
#' @rdname fishnet
#' @keywords bioregion
#' @importFrom terra rast crs<- as.polygons crs ext
#' @return A spatial vector polygon object of equal area grid cells covering
#' the defined area.
#' @references
#' Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum entropy
#' modeling of species geographic distributions. \emph{Ecological Modelling}
#' \strong{190}: 231-259.
#' @examples
#' d <- terra::vect(system.file("ex/nigeria.json", package="phyloregion"))
#' f <- fishnet(d, res = 0.75)
#' @export
fishnet <- function(mask, res = 0.5){
    s <- rast(res = res, ext(mask))
    crs(s) <- "epsg:4326"
    m <- as.polygons(s, dissolve = FALSE)
    m$grids <- paste0("v", seq_len(nrow(m)))
    m <- m[, "grids"]
    m[mask, ]
}
