progress <- function(x, FUN, ...) {
    env <- environment()
    pb_Total <- length(x)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3,
                         width = getOption("width")/2L)
    # wrapper around FUN
    wrapper <- function(...){
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        setTxtProgressBar(get("pb", envir=env), curVal +1)
        FUN(...)
    }
    r <- lapply(x, wrapper, ...)
    close(pb)
    r
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
#' @importFrom utils txtProgressBar setTxtProgressBar
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
#' @importFrom terra values values<- ncell setValues buffer extract crs
#' @importFrom terra project crs<- rasterize
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar object.size
#' @param trace Trace the function; trace = 2 or higher will be more voluminous.
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
polys2comm <- function(dat, res = 0.25, pol.grids = NULL, trace = 1, ...) {
    dat <- .matchnms(dat)

    if (!is.null(pol.grids)) {
        m <- pol.grids[, grep("grids", names(pol.grids)), drop = FALSE]
        pj <- "+proj=longlat +datum=WGS84"
        m <- suppressWarnings(invisible(project(m, "epsg:4326")))
        crs(m) <- pj

        s <- base::split(dat, f = dat$species)
        r <- rast(res = res, ext(m))
        ras <- rasterize(m, r, "grids", touches = TRUE)

        a <- as.data.frame(ras)
        b <- as.numeric(rownames(a))

        if (object.size(dat) > 15000L && interactive() && trace > 0) {
            f <- progress(s, function(x) {
                rs <- rasterize(x, r, touches = TRUE)
                k <- which(values(rs)==1)
                j <- as.character(intersect(k, b))
                unique(as.character(a[j,]))
            })
        } else {
            f <- lapply(s, function(x) {
                rs <- rasterize(x, r, touches = TRUE)
                k <- which(values(rs)==1)
                j <- as.character(intersect(k, b))
                unique(as.character(a[j,]))
            })
        }
        l <- lengths(f)
        dx <- cbind(species = rep(names(l), l),
                    grids = unlist(f, use.names = FALSE))
        dx <- na.omit(dx)
        y <- long2sparse(dx)
        g <- rowSums(y)

        i <- numeric(length(m))
        names(i) <- m$grids
        i[names(g)] <- g
        values(m) <- cbind(values(m), richness=i)
        z <- m[m$richness > 0, ]
    } else {
        s <- base::split(dat, f = dat$species)
        r <- rast(res=res, ext(dat))
        if (object.size(dat) > 15000L && interactive() && trace > 0) {
            files <- progress(s, function(x) {
                ras <- rasterize(x, r)
                which(values(ras)==1)
            })
        } else {
            files <- lapply(s, function(x) {
                ras <- rasterize(x, r)
                which(values(ras)==1)
            })
        }
        l <- lengths(files)
        row_nam <- factor(paste0("v", unlist(files)))
        y <- sparseMatrix(as.integer(row_nam), rep(seq_along(files), l), x = 1L,
                          dimnames = list(levels(row_nam), names(files)))
        g <- rowSums(y)
        values(r) <- paste0("v", seq_len(ncell(r)))
        i <- match(as.data.frame(r)[[1]], names(g))
        z <- setValues(r, g[i])
    }
    return(list(comm_dat = y, map = z))
}

#' @rdname rast2comm
#' @importFrom terra project vect crs<- rasterize as.data.frame
#' @importFrom terra values values<- ncell setValues buffer extract crs
#' @importFrom stats complete.cases predict
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
#' pt <- points2comm(dat = m, res = 0.5, lon = "lon", lat = "lat",
#'             species = "taxon") # This generates a list of two objects
#' head(pt[[1]])
#' @export
points2comm <- function(dat, res = 0.25, pol.grids = NULL, trace = 1, ...) {
    dat <- .matchnames(dat)
    dat <- na.omit(dat)

    dat <- vect(dat)

    if (!is.null(pol.grids)) {
        m <- pol.grids[, grep("grids", names(pol.grids)), drop = FALSE]
        pj <- "+proj=longlat +datum=WGS84"
        m <- suppressWarnings(invisible(project(m, "epsg:4326")))
        crs(m) <- pj

        s <- base::split(dat, f = dat$species)
        r <- rast(res = res, ext(m))
        ras <- rasterize(m, r, "grids", touches = TRUE)

        a <- as.data.frame(ras)
        b <- as.numeric(rownames(a))

        if (object.size(dat) > 15000L && interactive() && trace > 0) {
            f <- progress(s, function(x) {
                rs <- rasterize(x, r, touches = TRUE)
                k <- which(values(rs)==1)
                j <- as.character(intersect(k, b))
                as.character(a[j,])
            })
        } else {
            f <- lapply(s, function(x) {
                rs <- rasterize(x, r, touches = TRUE)
                k <- which(values(rs)==1)
                j <- as.character(intersect(k, b))
                as.character(a[j,])
            })
        }
        l <- lengths(f)
        dx <- cbind(species = rep(names(l), l),
                    grids = unlist(f, use.names = FALSE))
        dx <- na.omit(dx)
        y <- long2sparse(dx)
        tmp <- data.frame(grids=row.names(y), abundance=rowSums(y),
                          richness=rowSums(y>0))
        z <- merge(m, tmp, by = "grids")
    } else {
        s <- base::split(dat, f = dat$species)
        r <- rast(res=res, ext(dat))
        if (object.size(dat) > 15000L && interactive() && trace > 0) {
            files <- progress(s, function(x) {
                ras <- rasterize(x, r)
                which(values(ras)==1)
            })
        } else {
            files <- lapply(s, function(x) {
                ras <- rasterize(x, r)
                which(values(ras)==1)
            })
        }
        l <- lengths(files)
        row_nam <- factor(paste0("v", unlist(files)))
        y <- sparseMatrix(as.integer(row_nam), rep(seq_along(files), l), x = 1L,
                          dimnames = list(levels(row_nam), names(files)))
        g <- rowSums(y)
        values(r) <- paste0("v", seq_len(ncell(r)))
        i <- match(as.data.frame(r)[[1]], names(g))
        z <- setValues(r, g[i])
    }
    return(list(comm_dat = y, map = z))
}


