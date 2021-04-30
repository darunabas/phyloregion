make_poly <- function(file){
    if (!inherits(file, "RasterLayer")) file <- raster(file)
    pol <- rasterToPolygons(file, fun=NULL, dissolve=FALSE, na.rm=TRUE)
    suppressWarnings(invisible(proj4string(pol) <- crs(file)))
    pol$grids <- paste0("v", seq_len(nrow(pol)))
    xx <- as.data.frame(file, xy=TRUE, na.rm=TRUE)
    #Make dataframe of all xy coordinates
    xx$grids <- paste0("v", seq_len(nrow(xx)))
    m <- merge(pol, xx, by = "grids")
    m <- m[, c("grids", "x", "y")]
    m
}

foo <- function(file, rast=NULL) {
    if (!inherits(file, "RasterLayer")) file <- raster(file)
    if(!is.null(rast)) {
        if(!raster::compareRaster(file, rast)) stop("Raster objects are different")
    }
    y <- raster::as.data.frame(file, xy=TRUE, na.rm=TRUE, long=TRUE)
    y$grids <- paste0("v", seq_len(nrow(y)))
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
    r <- raster(ncol = 180, nrow = 180, resolution = res)
    extent(r) <- extent(p)
    r1 <- setValues(r, sample(x = 0:1, size = ncell(r), replace = TRUE))
    r1[!is.na(r1)] <- 0
    rp <- rasterize(x, r, field=1)
    res <- merge(rp, r1)
    names(res) <- x[[1]]
    return(res)
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
#' @param species a character string. The column with the species or taxon name.
#' Default = \dQuote{species}.
#' @param lon character with the column name of the longitude.
#' @param lat character with the column name of the latitude.
#' @param mask Only applicable to \code{points2comm}. If supplied, a polygon
#' shapefile covering the boundary of the survey region.
#' @param shp.grids if specified, the polygon shapefile of grid cells
#' with a column labeled \dQuote{grids}.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname raster2comm
#' @seealso \code{\link[mapproj]{mapproject}} for conversion of
#' latitude and longitude into projected coordinates system.
#' \code{\link{long2sparse}} for conversion of community data.
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values crs as.data.frame compareRaster
#' @importFrom raster rasterize setValues extent merge
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
    pol <- make_poly(files[1])
    pol <- pol[, "grids"]
    m <- progress(files, foo, rast=raster(files[1]))
    res <- do.call("rbind", m)
    y <- long2sparse(res)
    tmp <- data.frame(grids=row.names(y), richness=rowSums(y>0))
    z <- sp::merge(pol, tmp, by = "grids")
    z <- z[!is.na(z@data$richness), ]
    return(list(comm_dat = y, poly_shp = z))
}


#' @rdname raster2comm
#' @importFrom sp coordinates over CRS proj4string merge split merge
#' @importFrom sp spTransform
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar object.size
#' @importFrom raster raster res rasterize xyFromCell getValues
#' @examples
#' \donttest{
#' s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
#' sp <- random_species(100, species=5, shp=s)
#' pol <- polys2comm(dat = sp, species = "species")
#' head(pol[[1]])
#' }
#'
#' @export
polys2comm <- function(dat, res = 1, species = "species", shp.grids = NULL, ...) {

    dat <- dat[, species, drop = FALSE]
    names(dat) <- "species"

    if (!is.null(shp.grids)) {
        shp.grids <- shp.grids[, grepl("grids", names(shp.grids)), drop=FALSE]
        m <- shp.grids

        pj <- suppressWarnings(invisible(proj4string(m)[[1]]))
        suppressWarnings(invisible(proj4string(m) <- CRS(pj)))
        suppressWarnings(invisible(proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")))
        dat <- spTransform(dat, CRS(pj))

        spo <- sp::over(dat, m, returnList = TRUE)
        spo <- lapply(spo, unlist)
        ll <- lengths(spo)
        y <- data.frame(grids=unlist(spo), species=rep(dat@data$species, ll))

        y <- long2sparse(unique(y[, c("grids", "species")]))
        tmp <- data.frame(grids=row.names(y), richness=rowSums(y>0))
        z <- sp::merge(m, tmp, by = "grids")
        z <- z[!is.na(z@data$richness), ]
    } else {
        s <- split(dat, f = dat$species)
        files <- lapply(s, function(x) blank(x, res = res))
        pol <- make_poly(files[[1]])
        pol <- pol[, "grids"]
        m <- progress(files, foo, rast=raster(files[[1]]))
        spo <- do.call("rbind", m)
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
#' @importFrom raster extent
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
points2comm <- function(dat, mask = NULL, res = 1, lon = "decimallongitude",
                        lat = "decimallatitude", species = "species",
                        shp.grids = NULL, ...) {
    dat <- as.data.frame(dat)
    dat <- dat[, c(species, lon, lat)]
    names(dat) <- c("species", "lon", "lat")

    dat <- dat[complete.cases(dat), ]

    coordinates(dat) <- ~ lon + lat

    if (length(shp.grids) == 0) {
        if (is.null(mask)) {
            e <- raster::extent(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90))
            p <- as(e, "SpatialPolygons")
            m <- sp::SpatialPolygonsDataFrame(p, data.frame(sp = "x"))
            m <- fishnet(mask = m, res = res)
        } else (m <- fishnet(mask = mask, res = res))
    } else m <- shp.grids
    proj4string(dat) <- proj4string(m)
    x <- over(dat, m)
    y <- cbind(as.data.frame(dat), x)
    y <- y[complete.cases(y), ]
    Y <- long2sparse(y)
    tmp <- data.frame(grids=row.names(Y), abundance=rowSums(Y),
                      richness=rowSums(Y>0))
    z <- sp::merge(m, tmp, by = "grids")
    z <- z[!is.na(z@data$richness), ]
    return(list(comm_dat = Y, poly_shp = z))
}

