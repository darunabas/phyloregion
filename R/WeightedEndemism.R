#' Measure the distribution of narrow-ranged or endemic species.
#'
#' \code{weighted_endemism} is species richness inversely weighted
#' by species ranges.
#
#' @param x A (sparse) community matrix.
#' @keywords bioregion
#' @importFrom Matrix rowSums Diagonal Matrix colSums
#'
#' @return A data frame of species traits by site.
#'
#' @references
#' Crisp, M.D., Laffan, S., Linder, H.P. & Monro, A. (2001) Endemism in the
#' Australian flora. \emph{Journal of Biogeography} \strong{28}: 183–198.
#'
#' Daru, B.H., Farooq, H., Antonelli, A. & Faurby, S. (2020) Endemism
#' patterns are scale dependent. \emph{Nature Communications} \strong{11}
#' : 2115.
#'
#' @examples
#' library(raster)
#' data(africa)
#' Endm <- weighted_endemism(africa$comm)
#' m <- merge(africa$polys, data.frame(grids=names(Endm), WE=Endm), by="grids")
#' m <- m[!is.na(m@data$WE),]
#'
#' plot_swatch(m, values = m$WE,
#'             col = hcl.colors(20, palette = "Blue-Red 3", rev=FALSE))
#'
#' @export
weighted_endemism <- function(x){
  if(inherits(x, "matrix") && ncol(x)>2) x <- Matrix(x, sparse=TRUE)
  if(!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  x@x[x@x > 1e-8] <- 1  # we want to count species and not occurrences
  y <- rowSums(x %*% Diagonal(x = 1 / colSums(x) ) )
  y
}

