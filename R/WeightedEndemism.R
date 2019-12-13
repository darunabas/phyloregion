.SpecRich <- function(x){
  mm1 <- data.frame(table(x$grids))
  names(mm1) <- c("grids", "SR")
  mm1
}



#' Measure the distribution of narrow-ranged or endemic species.
#'
#' \code{weighted.endemism} is species richness inversely weighted
#' by species ranges.
#'
#' @param x A community matrix or data frame.
#' @param species If \code{format =} \dQuote{long} (the default),
#' the column with the species name.
#' @param grids The column with the sites or grids if \code{format =} \dQuote{long}.
#' @param format Format of the community composition data: \dQuote{long} or \dQuote{wide}
#' with species as columns and sites as rows.
#' @rdname weighted.endemism
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values
#' @importFrom sp CRS proj4string
#' @importFrom data.table as.data.table
#'
#' @return A data frame of species traits by site.
#'
#' @references
#' Crisp, M.D., Laffan, S., Linder, H.P. & Monro, A. (2001) Endemism in the Australian flora.
#' \emph{Journal of Biogeography} \strong{28}: 183–198.
#'
#' Laffan, S.W., & Crisp, M.D. (2003) Assessing endemism at multiple spatial scales,
#' with an example from the Australian vascular flora. \emph{Journal of Biogeography} \strong{30}: 511–520.
#' @examples
#' require(raster)
#' data(africa)
#'
#' Endm <- weighted.endemism(africa$comm)
#' m <- merge(africa$polys, Endm, by="grids")
#' m <- m[!is.na(m@data$WE),]
#'
#' plot_swatch(m, values = m$WE, k=20)
#'
#' @export weighted.endemism
weighted.endemism <- function(x, species="species", grids="grids", format="long"){

    if (format == "wide") {
        x <- data.frame(expand.grid(dimnames(provideDimnames(x)))[1:2],
                        as.vector(as.matrix(x)))
        x <- x[(x[, 3] > 0) & !is.na(x[, 3]), ]
        x <- x[sort.list(x[, 1]), ]
        x <- data.frame(grids = x[, 1], species = x[, 2])
        x <- as.data.frame(x)
        WE <- grids <- V1 <- NULL
        tmp <- .SpecRich(x)
        index <- match(x$grids, tmp$grids)
        SR <- tmp$SR[index]
        ff <- table(x$species)
        x$WE <- as.numeric(SR/ff[x$species]) # added as.number
        tmp <- as.data.table(x)
        res <- tmp[, sum(WE), by=grids]
        res <- res[, list(grids, WE=V1)]}
    else if (format == "long") {
        x <- as.data.frame(x)
        x <- x[, c(grids, species)]
        names(x) <- c("grids", "species")

        WE <- grids <- V1 <- NULL
        tmp <- .SpecRich(x)
        index <- match(x$grids, tmp$grids)
        SR <- tmp$SR[index]
        ff <- table(x$species)
        x$WE <- as.numeric(SR/ff[x$species]) # added as.number
        tmp <- as.data.table(x)
        res <- tmp[, sum(WE), by=grids]
        res <- res[, list(grids, WE=V1)]
    }
    res
}




