# write functions for each model
# model=c("RF", "GLM", "MAXENT", "GBM") to choose models
# run bioclim to generate points for species with limited sampling.

.arcp <- function (xy)
{
    if (nrow(xy) < 3)
        return(0)
    x.segmat <- cbind(xy, rbind(xy[2:nrow(xy), ], xy[1, ]))
    abs(sum(x.segmat[, 1] * x.segmat[, 4] - x.segmat[, 3] *
                x.segmat[, 2]))/2
}


.arcpspdf <- function (spdf)
{
    lar <- unlist(lapply(polygons(spdf)@polygons, function(x) unlist(lapply(x@Polygons,
                                                                            function(y) .arcp(y@coords)))))
    lhol <- unlist(lapply(polygons(spdf)@polygons, function(x) unlist(lapply(x@Polygons,
                                                                             function(y) y@hole))))
    sum(lar[!lhol]) - sum(lar[lhol])
}

MCP <- function (xy, percent = 95, unin = c("m", "km"), unout = c("ha",
                                                                  "km2", "m2"))
{
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should be of class SpatialPoints")
    if (ncol(coordinates(xy)) > 2)
        stop("xy should be defined in two dimensions")
    pfs <- proj4string(xy)
    if (length(percent) > 1)
        stop("only one value is required for percent")
    if (percent > 100) {
        warning("The MCP is estimated using all relocations (percent>100)")
        percent <- 100
    }
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy) != 1) {
            warning("xy should contain only one column (the id of the animals), id ignored")
            id <- factor(rep("a", nrow(as.data.frame(xy))))
        }
        else {
            id <- xy[[1]]
        }
    }
    else {
        id <- factor(rep("a", nrow(as.data.frame(xy))))
    }
    if (percent > 100) {
        warning("The MCP is estimated using all relocations (percent>100)")
        percent <- 100
    }
    if (min(table(id)) < 5)
        stop("At least 5 locations are required to fit a home range")
    id <- factor(id)
    xy <- as.data.frame(coordinates(xy))
    r <- split(xy, id)
    est.cdg <- function(xy) apply(xy, 2, mean)
    cdg <- lapply(r, est.cdg)
    levid <- levels(id)
    res <- SpatialPolygons(lapply(1:length(r), function(i) {
        k <- levid[i]
        df.t <- r[[levid[i]]]
        cdg.t <- cdg[[levid[i]]]
        dist.cdg <- function(xyt) {
            d <- sqrt(((xyt[1] - cdg.t[1])^2) + ((xyt[2] - cdg.t[2])^2))
            return(d)
        }
        di <- apply(df.t, 1, dist.cdg)
        key <- c(1:length(di))
        acons <- key[di <= quantile(di, percent/100)]
        xy.t <- df.t[acons, ]
        coords.t <- chull(xy.t[, 1], xy.t[, 2])
        xy.bord <- xy.t[coords.t, ]
        xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
        so <- Polygons(list(Polygon(as.matrix(xy.bord))), k)
        return(so)
    }))
    are <- unlist(lapply(1:length(res), function(i) {
        .arcpspdf(res[i, ])
    }))
    if (unin == "m") {
        if (unout == "ha")
            are <- are/10000
        if (unout == "km2")
            are <- are/1e+06
    }
    if (unin == "km") {
        if (unout == "ha")
            are <- are * 100
        if (unout == "m2")
            are <- are * 1e+06
    }
    df <- data.frame(id = unlist(lapply(1:nlevels(id), function(i) res[i]@polygons[[1]]@ID)),
                     area = are)
    row.names(df) <- df[, 1]
    res <- SpatialPolygonsDataFrame(res, df)
    if (!is.na(pfs))
        proj4string(res) <- CRS(pfs)
    return(res)
}


sampleBuffer <- function(x, n_points, width=2, limits=NULL){
    cc <- gBuffer(x, width=width)
    if(!is.null(limits)) cc <- crop(cc, limits)
    res <- spsample(cc, n_points, type="random")
    res
}


#' Species distribution models for a range of algorithms
#'
#' This function computes species distribution models using
#' four modelling algorithms: generalized linear models,
#' generalized boosted models, random forests, and maximum entropy (only if
#' \code{rJava} is available). Note: this is an experimental function, and
#' may change in the future.
#'
#' @param x A dataframe containing the species occurrences
#' and geographic coordinates. Column 1 labeled as "species", column 2 "lon",
#' column 3 "lat".
#' @param predictors RasterStack of environmental descriptors on which
#' the models will be projected
#' @param blank A blank raster upon which the prediction layer is aggregated to.
#' @param res Desired resolution of the predicted potential species distribution
#' (if blank raster is not specified).
#' @param pol A polygon shapefile specifying the boundary to restrict the
#' prediction. If not specified, a minimum convex polygon is estimated using
#' the input data frame of species occurrences.
#' @param tc Integer. Tree complexity. Sets the complexity of individual trees
#' @param lr Learning rate. Sets the weight applied to individual trees
#' @param bf Bag fraction. Sets the proportion of observations used in selecting
#' variables
#' @param n.trees Number of initial trees to fit. Set at 50 by default
#' @param k Number of groups
#' @param step.size Number of trees to add at each cycle
#' @param herbarium.rm Logical, remove points within 50 km of herbaria.
#' @param n_points Minimum number of points required to successfully run
#' a species distribution model
#' @rdname sdm
#' @importFrom raster values extent res<- crop extract predict resample merge
#' @importFrom raster stack calc buffer mask setValues extent<- nlayers maxValue
#' @importFrom raster crs<-
#' @importFrom sp coordinates CRS proj4string spsample HexPoints2SpatialPolygons
#' @importFrom sp Polygon SpatialPolygons polygons
#' @importFrom dismo kfold randomPoints evaluate threshold maxent gbm.step
#' @importFrom randomForest randomForest
#' @importFrom stats glm median formula
#' @importFrom grDevices chull
#' @importFrom rgeos gBuffer
#' @return A list with the following objects:
#' \itemize{
#'   \item \code{ensemble_raster} The ensembled raster that predicts
#'   the potential species distribution.
#'   \item \code{ensemble_AUC} The median AUCs of models.
#'   \item \code{data} The dataframe that was used to implement the model.
#'   \item \code{indiv_models} Raster layers for the separate models that
#'   predict the potential species distribution.
#'   \item \code{single_AUCs} The AUCs for the seperate models.
#' }
#' @references
#' Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum entropy
#' modeling of species geographic distributions. \emph{Ecological Modelling}
#' \strong{190}: 231-259.
#' @examples
#' \donttest{
#' library(raster)
#' # get predictor variables
#' f <- list.files(path=paste(system.file(package="phyloregion"), '/ex', sep=''),
#'                      pattern='.tif', full.names=TRUE )
#' preds <- stack(f)
#' #plot(preds)
#' # get species occurrences
#' d <- read.csv(system.file("ex/Bombax.csv", package="phyloregion"))
#'
#' # fit ensemble model for four algorithms
#' mod <- sdm(d, predictors = preds)
#' }

#' @export
sdm <- function(x, pol = NULL, predictors = NULL, blank = NULL, res = 1, tc = 2,
                lr = 0.001, bf = 0.75, n.trees = 50, step.size = n.trees, k=5,
                herbarium.rm = TRUE, n_points = 30) {
    x <- x[, c("species", "lon", "lat")]
    names(x) <- c("species", "lon", "lat")
    name.sp <- unique(x$species)

    if (is.null(predictors)) {
        stop("you need to specify RasterStack of environmental variables")
    }

    x <- x[, c("lon", "lat")]
    x <- x[!is.na(x$lon) & !is.na(x$lat),]
    x$source <- "raw"
    coordinates(x) <- ~lon+lat
    proj4string(x) <- CRS("+proj=longlat +datum=WGS84")

    # Remove points within 50 km of herbaria
    if (herbarium.rm) {
        dx <- read.csv(system.file("ex/IHfeb20.csv", package="phyloregion"))
        coordinates(dx) <- ~lon+lat
        proj4string(dx) <- CRS("+proj=longlat +datum=WGS84")
        herb_pol <- gBuffer(dx, width=0.5)
        x <- x[is.na(over(x, herb_pol)),]
    }

    fam_pol <- pol

    if (!is.null(pol)) {
        x <- x[pol,]
    }

    if (nrow(x) < 5) {
        e <- as(extent(x), "SpatialPolygons")
        crs(e) <- "+proj=longlat +datum=WGS84"
        pol <- suppressWarnings(invisible(gBuffer(e, width=2)))
    } else {
        pol <- suppressWarnings(invisible(gBuffer(MCP(x), width=2)))
    }
    #pol <- gBuffer(x, width=1)
    pol <- SpatialPolygonsDataFrame(pol, data.frame(id=1:length(pol)),
                                    match.ID = FALSE)

    x1 <- x[pol,]
    gc()

    if (is.null(blank)) {
        e <- raster::extent(pol)*1.5
        p <- as(e, "SpatialPolygons")
        r <- raster(ncol = 180, nrow = 180, resolution = res)
        extent(r) <- extent(p)
        r1 <- setValues(r, sample(x = 0:1, size = ncell(r), replace = TRUE))
        blank <- resample(predictors[[1]], r1)
        # set all values to zero
        blank[!is.na(blank)] <- 0
    }


    if (length(x1) > 0) {
        if(length(x1) < n_points){
            vv <- sampleBuffer(x1, n_points - length(x1), limits=pol)
            vv <- data.frame(vv)
            vv$source <- "random"
            names(vv)[c(1,2)] <- c("lon", "lat")
            x1 <- as.data.frame(x1)
            x1 <- rbind(x1, vv)
        }

        # OCCURRENCE POINTS
        occ_csv <- as.data.frame(x1)
        occ <- occ_csv[, c("lon", "lat")]
        fold <- dismo::kfold(occ, k=k)
        test <- 3
        pres_train <- occ[fold != test, ]
        pres_test <- occ[fold == test, ]

        if (!is.null(fam_pol)) {
            predictors <- raster::crop(predictors, fam_pol)
        } else {
            predictors <- raster::crop(predictors, pol)
        }

        gc()

        # only run if the maxent.jar file is available, in the right folder
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

        set.seed(10)
        backg <- dismo::randomPoints(mask=predictors, n=nrow(occ)*3,
                                     ext=extent(pol), extf = 1.1, warn = 0, p = occ)

        colnames(backg) <- c('lon', 'lat')

        set.seed(0)
        group <- dismo::kfold(backg, k=k)
        backg_train <- backg[group != test, ]
        backg_test <- backg[group == test, ]

        train <- rbind(pres_train, backg_train)
        pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
        envtrain <- extract(predictors, train)
        envtrain <- na.omit(data.frame(cbind(pa=pb_train, envtrain)))

        testpres <- data.frame(extract(predictors, pres_test) )
        testbackg <- data.frame(extract(predictors, backg_test))

        Preds <- names(envtrain)[-1]
        y <- names(envtrain)[1]
        Formula <- formula(paste(y," ~ ", paste(Preds, collapse=" + "),'+',
                           paste("I(", Preds, "^2)", sep = "",
                                 collapse = " + ")))

        # 1. Random Forest
        rf1 <- suppressWarnings(invisible(randomForest::randomForest(Formula, data=envtrain)))
        erf <- dismo::evaluate(testpres, testbackg, rf1)
        px <- raster::predict(predictors, rf1, ext=extent(pol))
        tr <- threshold(erf, 'spec_sens')
        trf <- raster::crop((px > tr), pol)
        # use the mask function
        zrf <- resample(mask(trf, pol), blank, method = "ngb")
        RF <- merge(zrf, blank)
        # TSS
        cm_rf <- data.frame(erf@confusion)
        tpr_rf <- cm_rf$tp / (cm_rf$tp + cm_rf$fn)
        tnr_rf <- cm_rf$tn / (cm_rf$fp + cm_rf$tn)
        tss_rf <- median(tpr_rf + tnr_rf - 1)
        gc()

        # 2. GLM
        gm <- suppressWarnings(invisible(glm(Formula, family = "binomial", data=envtrain)))
        ge <- dismo::evaluate(testpres, testbackg, gm)
        pg <- predict(predictors, gm, ext=extent(pol))
        gtr <- threshold(ge, 'spec_sens')
        gt <- raster::crop((pg > gtr), pol)
        # use the mask function
        gt1 <- resample(mask(gt, pol), blank, method = "ngb")
        GLM <- merge(gt1, blank)
        # TSS
        cm_gl <- data.frame(ge@confusion)
        tpr_gl <- cm_gl$tp / (cm_gl$tp + cm_gl$fn)
        tnr_gl <- cm_gl$tn / (cm_gl$fp + cm_gl$tn)
        tss_glm <- median(tpr_gl + tnr_gl - 1)
        gc()

        # 3. Maxent
        maxent_available <- FALSE
        if (maxent()) {
            maxent_available <- TRUE
            xm <- maxent(predictors, pres_train)
            xe <- evaluate(pres_test, backg_test, xm, predictors)
            mx <- predict(predictors, xm, ext=extent(pol), progress='')
            mr <- threshold(xe, 'spec_sens')
            mt <- raster::crop((mx > mr), pol)
            # use the mask function
            mt1 <- resample(mask(mt, pol), blank, method = "ngb")
            MX <- merge(mt1, blank)
            # TSS
            cm_mx <- data.frame(xe@confusion)
            tpr_mx <- cm_mx$tp / (cm_mx$tp + cm_mx$fn)
            tnr_mx <- cm_mx$tn / (cm_mx$fp + cm_mx$tn)
            tss_mx <- median(tpr_mx + tnr_mx - 1)
        }

        gc()

        # 4. Generalised boosted models, GBMs
        gbm_mod <- suppressWarnings(invisible(dismo::gbm.step(data = envtrain, gbm.x = Preds,
                                                              gbm.y=y, family = "bernoulli",
                                                              tree.complexity = tc,
                                                              learning.rate = lr,
                                                              bag.fraction = bf,
                                                              verbose = TRUE, silent = TRUE,
                                                              plot.main = FALSE,
                                                              step.size = step.size,
                                                              n.trees = n.trees)))

        egb <- suppressMessages(invisible(evaluate(testpres, testbackg, gbm_mod)))
        pred_gb <- raster::predict(object = predictors, model = gbm_mod,
                                   n.trees = gbm_mod$gbm.call$best.trees,
                                   type = "response", na.rm = TRUE)
        gb <- gbm_mod$cv.statistics$cv.threshold
        gbf <- raster::crop((pred_gb > gb), pol)
        # use the mask function
        gbm1 <- resample(mask(gbf, pol), blank, method = "ngb")
        GBM <- merge(gbm1, blank)
        # TSS
        cm_gb <- data.frame(egb@confusion)
        tpr_gb <- cm_gb$tp / (cm_gb$tp + cm_gb$fn)
        tnr_gb <- cm_gb$tn / (cm_gb$fp + cm_gb$tn)
        tss_gbm <- median(tpr_gb + tnr_gb - 1)
        gc()

        if (maxent_available) {
            models <- stack(RF, GLM, MX, GBM)
            names(models) <- c("RF", "GLM", "MAXENT", "GBM")
            m <- models[[which(!maxValue(models)==0)]]
            if (nlayers(m) > 1) {
                m <- calc(m, median, forceapply=TRUE)
            }
            spo <- m==1

            aucs <- c(auc_RF=erf@auc, auc_GLM=ge@auc, auc_MAXENT=xe@auc,
                      auc_GBM=egb@auc)
            m_auc <- median(aucs)

            indiv_TSS <- c(TSS_RF=tss_rf, TSS_GLM=tss_glm, TSS_MAXENT=tss_mx,
                           TSS_GBM=tss_gbm)
            TSS <- median(indiv_TSS)

            occ_csv$species <- name.sp
            occ_csv$AUC <- m_auc
            occ_csv$TSS <- TSS
            occ_csv <- occ_csv[, c("species", "lon", "lat", "source", "AUC",
                                   "TSS")]
        } else {
            models <- stack(RF, GLM, GBM)
            names(models) <- c("RF", "GLM", "GBM")
            m <- models[[which(!maxValue(models)==0)]]
            if (nlayers(m) > 1) {
                m <- calc(m, median, forceapply=TRUE)
            }
            spo <- m==1

            aucs <- c(auc_RF=erf@auc, auc_GLM=ge@auc, auc_GBM=egb@auc)
            m_auc <- median(aucs)

            indiv_TSS <- c(TSS_RF=tss_rf, TSS_GLM=tss_glm, TSS_GBM=tss_gbm)
            TSS <- median(indiv_TSS)

            occ_csv$species <- name.sp
            occ_csv$AUC <- m_auc
            occ_csv$TSS <- TSS
            occ_csv <- occ_csv[, c("species", "lon", "lat", "source", "AUC",
                                   "TSS")]
        }

        return(list(ensemble_raster = spo, ensemble_AUC = m_auc,
                    ensemble_TSS=TSS, data=occ_csv,
                    indiv_models=models, indiv_AUCs=aucs, indiv_TSS=indiv_TSS))
        gc()
    }
}





