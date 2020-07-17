make_poly <- function(file){
    dd <- raster(file)
    pol <- rasterToPolygons(dd, fun=NULL, dissolve=FALSE, na.rm=FALSE)
    proj4string(pol) <- crs(dd)
    pol$grids <- paste0("v", seq_len(nrow(pol)))
    xx <- as.data.frame(xyFromCell(dd, cell=seq_len(ncell(dd))))
    #Make dataframe of all xy coordinates
    xx$grids <- paste0("v", seq_len(nrow(xx)))
    m <- merge(pol, xx, by = "grids")
    #names(m)[2] <- "richness"
    names(m)[3] <- "lon"
    names(m)[4] <- "lat"
    m <- m[, c("grids", "lon", "lat")]
    m
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
#' @importFrom raster values crs
#' @importFrom sp CRS proj4string<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return
#' \itemize{
#'   \item comm_dat: (sparse) community matrix
#'   \item poly_shp: shapefile of grid cells with the values per cell.
#' }
#' @examples
#' \donttest{
#' fdir <- system.file("NGAplants", package="phyloregion")
#' files <- file.path(fdir, dir(fdir))
#' ras <- raster2comm(files)
#' head(ras[[1]])
#' }
#'
#' @export
raster2comm <- function(files) {
    poly <- make_poly(files[1])
    tmp <- raster(files[1])
    fg <- as.data.frame(xyFromCell(tmp, cell=1:ncell(tmp)))
    ind1 <- paste(as.character(fg$x), as.character(fg$y), sep="_")
    ind2 <- paste(as.character(poly$lon), as.character(poly$lat), sep="_")
    index <- match(ind1, ind2)
    res <- NULL
    cells <- as.character(poly$grids)
    if (interactive()){
        pb <- txtProgressBar(min = 0, max = length(files), style = 3,
                             width = getOption("width")/2L)
    }

    for(i in seq_along(files)) {
        obj <- raster(files[i])
        tmp <- values(obj)
        tmp <- cells[index[tmp > 0]]
        tmp <- tmp[!is.na(tmp)]
        if(length(tmp) > 0) res <- rbind(res, cbind(tmp, names(obj)))
        if (interactive()) setTxtProgressBar(pb, i)
    }
    if(length(tmp) > 0) colnames(res) <- c("grids", "species")
    y <- long2sparse(as.data.frame(res))
    return(list(comm_dat = y, poly_shp = poly))
}


#' @rdname raster2comm
#' @importFrom sp coordinates over CRS proj4string merge split merge
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar object.size
#' @importFrom raster raster res rasterize xyFromCell getValues
#' @param trace Trace the function; trace = 2 or higher will be more voluminous.
#' @examples
#' \donttest{
#' s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
#' sp <- random_species(100, species=5, shp=s)
#' pol <- polys2comm(dat = sp, species = "species", trace=0)
#' head(pol[[1]])
#' }
#'
#' @export
polys2comm <- function(dat, res = 1, species = "species", trace = 1, ...) {

    dat <- dat[, species, drop = FALSE]
    names(dat) <- "species"

    e <- raster(dat)
    res(e) <- res
    s <- split(dat, f = dat$species)
    w <- rasterize(s[[1]], e)
    poly <- make_poly(w)
    fg <- as.data.frame(xyFromCell(w, cell = 1:ncell(w)))
    ind1 <- paste(as.character(fg$x), as.character(fg$y), sep = "_")
    ind2 <- paste(as.character(poly$lon), as.character(poly$lat), sep = "_")
    index <- match(ind1, ind2)
    r <- NULL
    cells <- as.character(poly$grids)[index]
    if (object.size(dat) > 150000L && interactive() && trace > 0) {
        m <- progress(s, function(x) {
            obj <- rasterize(x, e)
            tmp <- getValues(obj)
            cells[!is.na(tmp) & (tmp>0)]
        })
    } else {
        m <- lapply(s, function(x) {
            obj <- rasterize(x, e)
            tmp <- getValues(obj)
            cells[!is.na(tmp) & (tmp>0)]
        })
    }

    spo <- data.frame(grids = unlist(m), species = rep(labels(s), lengths(m)))
    y <- long2sparse(spo)

    z <- data.frame(grids = row.names(y), richness = rowSums(y > 0))
    z <- sp::merge(poly, z, by = "grids")
    z <- z[!is.na(z@data$richness), ]
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
#'             species = "taxon")
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
        if (length(mask) == 0) {
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

