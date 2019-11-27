sampleBuffer <- function(x, n_points=20, width=2, limits=NULL){
  cc <- gBuffer(x, width=width)
  if(!is.null(limits)) cc <- crop(cc, limits)
  res <- spsample(cc, n_points, type="random")
  res
}



#' Global scale species distribution modeling
#'
#' The \code{sdm} function enables the process in which species occurrence
#' data and ecological and environmental variables are modeled together to
#' predict species ranges.
#'
#' @param files A community matrix
#' @rdname sdm
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values
#' @importFrom sp CRS proj4string<-
#'
#' @export
#' @return
#' \item{psim}{A site × site phylogenetic beta diversity distance matrix}
#'
#' @references
#'
#' \insertRef{Philips2006}{bioregion}
#'
#' @examples
#' fdir <- system.file("Aloes", package="bioregion")
#' files <- file.path(fdir, dir(fdir))
#' res <- raster2comm(files)
sdm <- function(x) {
  name.sp <- S[i]
  newdat <- subset(mydata5, mydata5$species==name.sp)

  ###############################################################################
  dups2 <- duplicated(newdat[, c("species", "lon", "lat")])
  b <- newdat[!dups2,]

  # we only need columns 2 and 3:
  b <- b[,2:3]
  b <- b[!is.na(b$lon) & !is.na(b$lat),]
  coordinates(b)=~lon+lat
  proj4string(b) = CRS("+proj=longlat +datum=WGS84")
  b1 <- b[s,]

  if(length(b1) >0){
    if(length(b1) < 20){
      vv <- sampleBuffer(b1, 20 - length(b1), limits=s)
      b1 <- vv+b1
      b1 <- data.frame(b1)
      names(b1) <- c("lon", "lat")
    }

    # remove NA rows
    occ <- b1[!is.na(b1$lon) & !is.na(b1$lat),]
    ###<-----find ways of saving this as a point data.

    occ <- data.frame(occ)
    occ_csv <- occ
    occ_csv$name <- name.sp
    occ_csv <- data.frame(occ_csv)
    # only run if the maxent.jar file is available, in the right folder
    jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

    # OCCURRENCE POINTS
    fold <- kfold(occ, k=2)
    pres_train <- occ[fold != 1, ]
    pres_test <- occ[fold == 1, ]

    backg <- randomPoints(predictors, n=1000, ext=extent(s), extf = 1.25)
    colnames(backg) = c('lon', 'lat')
    group <- kfold(backg, k=2)
    backg_train <- backg[group != 1, ]
    backg_test <- backg[group == 1, ]

    xm <- maxent(predictors, pres_train)

    e1 <- evaluate(pres_test, backg_test, xm, predictors)
    px <- predict(predictors, xm, ext=extent(s), progress='')
    tr <- threshold(e1, 'spec_sens')
    ttt <- (px > tr)
    ttt1 <- crop(ttt, s)
    # now use the mask function
    ttt2 <- mask(ttt1, s)
    z <- resample(ttt2, dd, method = "ngb")
    w <- merge(z, dd)

    # LEVEL 2
    pol_1 <- rasterToPolygons(w, fun=function(x){x>0}, dissolve=TRUE)
    proj4string(pol_1) = CRS("+proj=longlat +datum=WGS84")

    occ_csv1 <- occ_csv
    coordinates(occ_csv1)=~lon+lat
    proj4string(occ_csv1) = CRS("+proj=longlat +datum=WGS84")
    cc <- gBuffer(occ_csv1, width=10)
    pol <- gIntersection(cc, pol_1, byid=TRUE)
    proj4string(pol) = CRS("+proj=longlat +datum=WGS84")

    tttt1 <- crop(ttt, pol)
    # now use the mask function
    tttt2 <- mask(tttt1, pol)
    zz <- resample(tttt2, dd, method = "ngb")
    ww <- merge(zz, dd)
    ww
    #writeRaster(ww, filename=paste("rasters/", name.sp, sep=""), format='GTiff', overwrite=TRUE)
    #write.table(occ_csv, paste0("CSVs/", name.sp, ".csv"), col.names = TRUE, quote = FALSE, row.names = F, sep = "\t")
    #write.csv(occ_csv, filename=paste0("CSVs/", name.sp, ".csv"), row.names = F)
  }
}


