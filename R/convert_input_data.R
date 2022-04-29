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
#' The functions \code{points2comm}, \code{polys2comm}, \code{raster2comm}
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
#' @rdname raster2comm
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
#' ras <- raster2comm(files) # Note, this function generates
#'      # a list of two objects
#' head(ras[[1]])
#' }
#'
#' @export
raster2comm <- function(files) {
    r <- raster::raster(files[1])
    m <- progress(files, foo, rast=raster::raster(files[1]))
    res <- do.call("rbind", m)
    if(!(nrow(res) > 0)) stop("Raster files probably empty!")
    y <- long2sparse(res)
    tmp <- data.frame(grids=row.names(y), richness=rowSums(y>0))
    if(raster::ncell(r) > 10000L) {
        r <- raster::raster(files[1])
        xy <- raster::as.data.frame(r, na.rm = FALSE, long = TRUE)
        xy$grids <- paste0("v", seq_len(nrow(xy)))
        xy <- na.omit(xy)
        xy <- xy[, "grids", drop = FALSE]
        xy$richness <- 0L
        ind1 <- rbind(tmp, xy)
        ind1 = ind1[!duplicated(ind1$grids), ]

        r[1:raster::ncell(r)] <- paste0("v", seq_len(raster::ncell(r)))
        index <- match(raster::values(r), ind1$grids)
        z <- raster::setValues(r, ind1$richness[index])
    } else {
        pol <- make_poly(files[1]) # no makepoly <<----
        pol <- pol[, "grids"]
        z <- sp::merge(pol, tmp, by = "grids")
        z <- z[!is.na(z@data$richness), ]
    }
    return(list(comm_dat = y, poly_shp = z))
}


#' @rdname raster2comm
#' @importFrom sp coordinates over CRS proj4string merge split merge
#' @importFrom sp spTransform
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar object.size
#' @param trace Trace the function; trace = 2 or higher will be more voluminous.
#' @examples
#' \donttest{
#' s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
#' sp <- random_species(100, species=5, shp=s)
#' pol <- polys2comm(dat = sp, species = "species")
#' head(pol[[1]])
#' }
#'
#' @export
polys2comm <- function(dat, res = 1, shp.grids = NULL, trace = 1,...) {

    dat <- .matchnms(dat)

    if (!is.null(shp.grids)) {
        shp.grids <- shp.grids[, grepl("grids", names(shp.grids)), drop=FALSE]
        m <- shp.grids

        pj <- suppressWarnings(invisible(proj4string(m)[[1]]))
        suppressWarnings(invisible(proj4string(m) <- CRS(pj)))
        suppressWarnings(invisible(proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")))
        dat <- spTransform(dat, CRS(pj))

        s <- split(dat, f = dat$species)

        if (object.size(dat) > 150000L && interactive() && trace > 0) {
            f <- progress(s, function(x) {
                tryCatch({
                spo <- sp::over(x, m, returnList = TRUE)
                spo <- lapply(spo, unlist)
                ll <- lengths(spo)
                data.frame(grids=unlist(spo), species=rep(x@data$species, ll))
                }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
            })
            result <- Filter(Negate(is.null), f)
            y <- do.call(rbind, result)

        } else {
            f <- lapply(s, function(x) {
                tryCatch({
                spo <- sp::over(x, m, returnList = TRUE)
                spo <- lapply(spo, unlist)
                ll <- lengths(spo)
                data.frame(grids=unlist(spo), species=rep(x@data$species, ll))
                }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
            })
            result <- Filter(Negate(is.null), f)
            y <- do.call(rbind, result)
        }

        y <- long2sparse(unique(y[, c("grids", "species")]))
        tmp <- data.frame(grids=row.names(y), richness=rowSums(y>0))
        z <- sp::merge(m, tmp, by = "grids")
        z <- z[!is.na(z@data$richness), ]
    } else {
        s <- split(dat, f = dat$species)
        if (object.size(dat) > 150000L && interactive() && trace > 0) {
            files <- progress(s, function(x) blank(x, res = res))
            pol <- make_poly(files[[1]])
            pol <- pol[, "grids"]
            m <- progress(files, foo, rast=raster::raster(files[[1]]))
            spo <- do.call("rbind", m)
        } else {
            files <- lapply(s, function(x) blank(x, res = res))
            pol <- make_poly(files[[1]])
            pol <- pol[, "grids"]
            m <- lapply(files, foo, rast=raster::raster(files[[1]]))
            spo <- do.call("rbind", m)
        }
        y <- long2sparse(spo)
        tmp <- data.frame(grids=row.names(y), richness=rowSums(y>0))
        z <- sp::merge(pol, tmp, by = "grids")
        z <- z[!is.na(z@data$richness), ]
    }

    return(list(comm_dat = y, poly_shp = z))
}


#' @rdname raster2comm
#' @importFrom sp coordinates<- over CRS proj4string merge
#' @importFrom sp coordinates<- CRS proj4string<- SpatialPolygonsDataFrame
#' @importFrom stats complete.cases
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

