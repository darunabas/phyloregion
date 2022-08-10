make_poly <- function(file){
    if (!inherits(file, "RasterLayer")) file <- raster::raster(file)
    pol <- raster::rasterToPolygons(file, fun=NULL, dissolve=FALSE, na.rm=FALSE)
    suppressWarnings(invisible(proj4string(pol) <- raster::crs(file)))
    pol$grids <- paste0("v", seq_len(nrow(pol)))
    xx <- raster::as.data.frame(file, xy=TRUE, na.rm=FALSE)
    #Make dataframe of all xy coordinates
    xx$grids <- paste0("v", seq_len(nrow(xx)))
    m <- raster::merge(pol, xx, by = "grids")
    m <- m[, c("grids", "x", "y")]
    m
}


foo <- function(file, rast=NULL) {
    if (!inherits(file, "RasterLayer")) file <- raster::raster(file)
    if(!is.null(rast)) {
        if(!raster::compareRaster(file, rast)) stop("Raster objects are different")
    }
    y <- raster::as.data.frame(file, xy=TRUE, na.rm=FALSE, long=TRUE)
    y$grids <- paste0("v", seq_len(nrow(y)))
    y <- na.omit(y)
    y <- y[y$value>0,, drop=FALSE]
    y <- y[, c("grids", "layer")]
    names(y)[2] <- "species"
    return(y)
}

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

blank <- function(x, res=NULL) {
    e <- raster::extent(c(-180, 180, -90, 90))
    p <- as(e, "SpatialPolygons")
    r <- raster::raster(ncol = 180, nrow = 180, resolution = res)
    raster::extent(r) <- raster::extent(p)
    r1 <- raster::setValues(r, sample(x = 0:1, size = raster::ncell(r), replace = TRUE))
    r1[!is.na(r1)] <- 0
    rp <- raster::rasterize(x, r, field=1)
    res <- raster::merge(rp, r1)
    names(res) <- x[[1]]
    return(res)
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
#' The functions \code{points2comm}, \code{vect2comm}, \code{rast2comm}
#' provide convenient interfaces to convert raw distribution data often
#' available as point records, polygons and raster layers,
#' respectively, to a community composition data frame at varying spatial grains
#' and extents for downstream analyses.
#'
#' @param files  list of raster layer objects with the same spatial
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
#' @param mask Only applicable to \code{points2comm}. If supplied, a polygon
#' shapefile covering the boundary of the survey region.
#' @param shp.grids if specified, the polygon shapefile of grid cells
#' with a column labeled \dQuote{grids}.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname rast2comm
#' @seealso \code{\link[mapproj]{mapproject}} for conversion of
#' latitude and longitude into projected coordinates system.
#' \code{\link{long2sparse}} for conversion of community data.
#' @importFrom sp CRS proj4string<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return Each of these functions generate a list of two objects as follows:
#' \itemize{
#'   \item comm_dat: (sparse) community matrix
#'   \item poly_shp: shapefile of grid cells with the values per cell.
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
#' @importFrom sf st_intersects
#' @importFrom sp coordinates over CRS proj4string merge split merge
#' @importFrom sp spTransform
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar object.size
#' @param trace Trace the function; trace = 2 or higher will be more voluminous.
#' @examples
#' \donttest{
#' s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
#' sp <- random_species(100, species=5, shp=s)
#' pol <- vect2comm(dat = sp, species = "species")
#' head(pol[[1]])
#' }
#'
#' @export
vect2comm <- function(dat, res = 1, shp.grids = NULL, trace = 1,...) {
    dat <- .matchnms(dat)

    if (!is.null(shp.grids)) {
        shp.grids <- shp.grids[, grep("grids", names(shp.grids)), drop=FALSE]
        m <- shp.grids
        pj <- crs(m, proj=TRUE)
        m <- project(m, pj)
        crs(dat) <- "+proj=longlat +datum=WGS84"
        dat <- project(dat, pj)
        dat <- sf::st_as_sf(dat)
        m <- sf::st_as_sf(m)
        s <- base::split(dat, f = dat$species)

        if (object.size(dat) > 150000L && interactive() && trace > 0) {
            f <- progress(s, function(x) {
                tryCatch({
                    l <- as.data.frame(sf::st_join(x, m, join = st_intersects))
                    l[-3]
                }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
            })
            result <- Filter(Negate(is.null), f)
            y <- do.call(rbind, result)

        } else {
            f <- lapply(s, function(x) {
                tryCatch({
                    l <- as.data.frame(sf::st_join(x, m, join = st_intersects))
                    l[-3]
                }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
            })
            result <- Filter(Negate(is.null), f)
            y <- do.call(rbind, result)
        }

        tmp <- long2sparse(y)
        g <- rowSums(tmp)
        i <- match(m$grids, names(g))
        tmp <- cbind(m, richness=g[i])
        z <- vect(tmp) # weird addition!

    } else {
        s <- base::split(dat, f = dat$species)
        r <- rast(res=res, ext(dat))
        if (object.size(dat) > 150000L && interactive() && trace > 0) {
            files <- progress(s, function(x) {
                ras <- rasterize(x, r)
                which(values(ras)==1)
            })
        } else {
            files <- progress(s, function(x) {
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

    return(list(comm_dat = y, raster = z))
}

#' @rdname rast2comm
#' @importFrom sp coordinates<- over CRS proj4string merge
#' @importFrom sp coordinates<- CRS proj4string<- SpatialPolygonsDataFrame
#' @importFrom stats complete.cases predict
#' @examples
#' s <- readRDS(system.file("nigeria/nigeria.rds", package = "phyloregion"))
#'
#' set.seed(1)
#' m <- data.frame(sp::spsample(s, 10000, type = "nonaligned"))
#' names(m) <- c("lon", "lat")
#' species <- paste0("sp", sample(1:1000))
#' m$taxon <- sample(species, size = nrow(m), replace = TRUE)
#'
#' pt <- points2comm(dat = m, mask = s, res = 0.5, lon = "lon", lat = "lat",
#'             species = "taxon") # Note, this generates a list of two objects
#' head(pt[[1]])
#' @export
points2comm <- function(dat, mask = NULL, res = 1, shp.grids = NULL, ...) {
    dat <- .matchnames(dat)

    #dat <- raster::as.data.frame(dat)
    #dat <- dat[, c(species, lon, lat)]
    #names(dat) <- c("species", "lon", "lat")

    dat <- na.omit(dat)

    coordinates(dat) <- ~ lon + lat

    if (length(shp.grids) == 0) {
        if (is.null(mask)) {
            e <- raster::extent(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90))
            p <- as(e, "SpatialPolygons")
            m <- sp::SpatialPolygonsDataFrame(p, data.frame(sp = "x"))
            m <- fishnet(mask = m, res = res)
            proj4string(dat) <- proj4string(m)
        } else {
            m <- fishnet(mask = mask, res = res)
            proj4string(dat) <- proj4string(m)
        }
    } else {
        shp.grids <- shp.grids[, grepl("grids", names(shp.grids)), drop=FALSE]
        m <- shp.grids
        pj <- suppressWarnings(invisible(proj4string(m)[[1]]))
        suppressWarnings(invisible(proj4string(m) <- CRS(pj)))
        suppressWarnings(invisible(proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")))
        dat <- spTransform(dat, CRS(pj))
    }

    #proj4string(dat) <- proj4string(m)
    x <- over(dat, m)
    y <- cbind(raster::as.data.frame(dat), x)
    y <- na.omit(y)
    Y <- long2sparse(y)
    tmp <- data.frame(grids=row.names(Y), abundance=rowSums(Y),
                      richness=rowSums(Y>0))
    z <- sp::merge(m, tmp, by = "grids")
    z <- z[!is.na(z@data$richness), ]
    return(list(comm_dat = Y, poly_shp = z))
}

