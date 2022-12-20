#progress <- function(x, FUN, ...) {
#    env <- environment()
#    pb_Total <- length(x)
#    counter <- 0
#    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3,
#                         width = getOption("width")/2L)
    # wrapper around FUN
#    wrapper <- function(...){
#        curVal <- get("counter", envir = env)
#        assign("counter", curVal +1 ,envir=env)
#        setTxtProgressBar(get("pb", envir=env), curVal +1)
#        FUN(...)
#    }
#    r <- lapply(x, wrapper, ...)
#    close(pb)
#    r
#}


.matchnames <- function(x) {
    x <- as.data.frame(x)
    nat <- colnames(x)
    X <- paste(c("\\blongitude\\b", "\\bdecimalLongitude\\b", "\\bLong\\b",
                 "\\bx\\b", "\\blon\\b"), collapse = "|")
    Y <- paste(c("\\blatitude\\b", "\\bdecimalLatitude\\b", "\\bLati\\b",
                 "\\by\\b", "\\lat\\b"), collapse = "|")
    SP <- paste(c("\\bspecies\\b", "\\bbinomial\\b", "\\bbinomil\\b",
                  "\\btaxon\\b"), collapse = "|")

    lon <- nat[grepl(X, nat, ignore.case = TRUE)]
    lat <- nat[grepl(Y, nat, ignore.case = TRUE)]
    species <- nat[grepl(SP, nat, ignore.case = TRUE)]
    x <- x[, c(species, lon, lat)]
    names(x) <- c("species", "lon", "lat")
    return(x)
}


.matchnms <- function(x) {
    nat <- names(x)
    SP <- paste(c("\\bspecies\\b", "\\bbinomial\\b", "\\bbinomil\\b",
                  "\\btaxon\\b"), collapse = "|")
    species <- nat[grepl(SP, nat, ignore.case = TRUE)]
    x <- x[, species, drop = FALSE]
    names(x) <- "species"
    return(x)
}

#' Convert raw input distribution data to community
#'
#' The functions \code{points2comm}, \code{polys2comm}, \code{rast2comm}
#' provide convenient interfaces to convert raw distribution data often
#' available as point records, polygons and raster layers,
#' respectively, to a community composition data frame at varying spatial grains
#' and extents for downstream analyses.
#'
#' @param files  list of SpatRaster layer objects with the same spatial
#' extent and resolution.
#' @param dat layers of merged maps corresponding to species polygons for
#' \code{polys2comm}; or point occurrence data frame for \code{points2comm},
#' with at least three columns:
#' \itemize{
#'   \item Column 1: \code{species} (listing the taxon names)
#'   \item Column 2: \code{decimallongitude} (corresponding to decimal longitude)
#'   \item Column 3: \code{decimallatitude} (corresponding to decimal latitude)
#' }
#' @param res the grain size of the grid cells in decimal degrees (default).
#' @param pol.grids if specified, the vector polygon of grid cells
#' with a column labeled \dQuote{grids}.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname rast2comm
#' @seealso \code{\link[mapproj]{mapproject}} for conversion of
#' latitude and longitude into projected coordinates system.
#' \code{\link{long2sparse}} for conversion of community data.
#' @return Each of these functions generate a list of two objects as follows:
#' \itemize{
#'   \item comm_dat: (sparse) community matrix
#'   \item map: vector or raster of grid cells with the values per cell for
#'   mapping.
#' }
#' @examples
#' \donttest{
#' fdir <- system.file("NGAplants", package="phyloregion")
#' files <- file.path(fdir, dir(fdir))
#' ras <- rast2comm(files) # Note, this function generates
#'      # a list of two objects
#' head(ras[[1]])
#' }
#'
#' @export
rast2comm <- function(files) {
    x <- terra::rast(files)
    ind <- which(values(x)==1, arr.ind = TRUE)
    row_nam <- factor(paste0("v", ind[,1]))
    y <- sparseMatrix(as.integer(row_nam), ind[, 2], x = 1L,
                      dimnames = list(levels(row_nam), names(x)))
    blank <- x[[1]]
    g <- rowSums(y)
    values(blank) <- paste0("v", seq_len(ncell(blank)))
    i <- match(as.data.frame(blank)[[1]], names(g))
    z <- setValues(blank, g[i])
    return(list(comm_dat = y, raster = z))
}

#' @rdname rast2comm
#' @importFrom terra values values<- ncell setValues buffer crs
#' @importFrom terra project crs<- rasterize ext
#' @importFrom methods as
#' @examples
#' \donttest{
#' require(terra)
#' s <- vect(system.file("ex/nigeria.json", package="phyloregion"))
#' sp <- random_species(100, species=5, pol=s)
#' pol <- polys2comm(dat = sp)
#' head(pol[[1]])
#' }
#'
#' @export
polys2comm <- function(dat, res = 0.25, pol.grids = NULL, ...) {
    dat <- .matchnms(dat)

    if (!is.null(pol.grids)) {
        m <- pol.grids[, grep("grids", names(pol.grids)), drop = FALSE]
    } else {
        m <- fishnet(mask = ext(dat), res = res)
        crs(m) <- "+proj=longlat +datum=WGS84"
    }
    pj <- "+proj=longlat +datum=WGS84"
    m <- suppressWarnings(invisible(project(m, "epsg:4326")))
    crs(m) <- pj
    j <- relate(dat, m, "intersects", pairs=TRUE)
    y <- m$grids[j[,2]]
    spp <- dat$species[j[,1]]
    r <- data.frame(species=spp, grids=y)
    y <- long2sparse(r)
    tmp <- data.frame(grids=row.names(y), abundance=rowSums(y),
                      richness=rowSums(y>0))
    z <- merge(m, tmp, by = "grids")
    return(list(comm_dat = y, map = z))
}

#' @rdname rast2comm
#' @importFrom terra project vect crs<- rasterize as.data.frame intersect
#' @importFrom terra values values<- ncell setValues buffer crs merge
#' @importFrom terra relate ext
#' @importFrom stats  predict
#' @importFrom methods as
#' @examples
#' library(terra)
#' s <- vect(system.file("ex/nigeria.json", package="phyloregion"))
#' set.seed(1)
#' m <- as.data.frame(spatSample(s, 1000, method = "random"),
#'                    geom = "XY")[-1]
#' names(m) <- c("lon", "lat")
#' species <- paste0("sp", sample(1:100))
#' m$taxon <- sample(species, size = nrow(m), replace = TRUE)
#'
#' pt <- points2comm(dat = m, res = 0.5) # This generates a list of two objects
#' head(pt[[1]])
#' @export
points2comm <- function(dat, res = 0.25, pol.grids = NULL, ...) {
    dat <- .matchnames(dat)
    dat <- na.omit(dat)

    if (!is.null(pol.grids)) {
        m <- pol.grids[, grep("grids", names(pol.grids)), drop = FALSE]
    } else {
        m <- fishnet(mask = ext(vect(dat)), res = res)
        crs(m) <- "+proj=longlat +datum=WGS84"
    }
    pj <- "+proj=longlat +datum=WGS84"
    m <- suppressWarnings(invisible(project(m, "epsg:4326")))
    crs(m) <- pj
    b <- vect(dat)
    j <- relate(b, m, "intersects", pairs=TRUE)
    y <- m$grids[j[,2]]
    spp <- b$species[j[,1]]
    r <- data.frame(species=spp, grids=y)
    y <- long2sparse(r)
    tmp <- data.frame(grids=row.names(y), abundance=rowSums(y),
                      richness=rowSums(y>0))
    z <- merge(m, tmp, by = "grids")
    return(list(comm_dat = y, map = z))
}
