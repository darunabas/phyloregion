#' Measure the distribution of narrow-ranged or endemic species.
#'
#' \code{weighted_endemism} is species richness inversely weighted
#' by species ranges.
#
#' @param x A (sparse) community matrix.
#' @rdname weighted.endemism
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values
#' @importFrom sp CRS proj4string
#' @importFrom Matrix rowSums
#'
#' @return A data frame of species traits by site.
#'
#' @references
#' Crisp, M.D., Laffan, S., Linder, H.P. & Monro, A. (2001) Endemism in the Australian flora.
#' \emph{Journal of Biogeography} \strong{28}: 183–198.
## Laffan, S.W., & Crisp, M.D. (2003) Assessing endemism at multiple spatial scales,
## with an example from the Australian vascular flora. \emph{Journal of Biogeography} \strong{30}: 511–520.
#' @examples
#' require(raster)
#' data(africa)
#' Endm <- weighted_endemism(africa$comm)
#' m <- merge(africa$polys, data.frame(grids=names(Endm), WE=Endm), by="grids")
#' m <- m[!is.na(m@data$WE),]
#'
#' plot_swatch(m, values = m$WE, k=20)
#'
#' @export
weighted_endemism <- function(x){
  if(inherits(x, "matrix") && ncol(x)>2) x <- Matrix(x, sparse=TRUE)
  if(!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  x@x[x@x > 1] <- 1  # we want to count species and not occurrences
  rowSums(x %*% Diagonal(x = 1 / colSums(x) ) )
}

