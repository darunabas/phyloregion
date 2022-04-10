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

# write functions for each model
# model=c("RF", "GLM", "MAXENT", "GBM") to choose models
# run bioclim to generate points for species with limited sampling.


#sampleBuffer <- function(x, n_points, width=2, limits=NULL){
#    cc <- gBuffer(x, width=width)
#    if(!is.null(limits)) cc <- crop(cc, limits)
#    res <- spsample(cc, n_points, type="random")
#    res
#}

#.blank_raster <- function(res) {
#    e <- raster::extent(c(-180, 180, -90, 90))
#    p <- as(e, "SpatialPolygons")
#    r <- raster(ncol = 180, nrow = 180, resolution = res)
#    extent(r) <- extent(p)
#    blank <- setValues(r, sample(x = 0:1, size = ncell(r), replace = TRUE))
#    # set all values to zero
#    blank[!is.na(blank)] <- 0
#    return(blank)
#}

.more_points <- function(pts, preds) {
    e <- extent(pts)
    x <- as.data.frame(pts)
    bc <- dismo::bioclim(preds, pts)
    p <- dismo::predict(preds, bc, tail='high', ext = e)
    set.seed(20220410)
    vv <- suppressWarnings(invisible(as.data.frame(randomPoints(mask=p,
                            n=nrow(x), prob=TRUE, ext = e))))
    vv$source <- "random"
    names(vv)[c(1,2)] <- c("lon", "lat")
    res <- rbind(x, vv)
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
#' @param n.points Minimum number of points required to successfully run
#' a species distribution model
#' @param prob_method Transform probabilities of presence into presence and
#' absence ("presab", the default) for binary transformation or "raw" for
#' overall projection.
#' @rdname sdm
#' @importFrom raster values extent res<- crop extract predict resample merge
#' @importFrom raster stack calc buffer mask setValues extent<- nlayers maxValue
#' @importFrom raster crs<-
#' @importFrom sp coordinates CRS proj4string spsample HexPoints2SpatialPolygons
#' @importFrom sp Polygon SpatialPolygons polygons
#' @importFrom dismo kfold randomPoints evaluate threshold maxent gbm.step
#' @importFrom randomForest randomForest
#' @importFrom stats glm median formula gaussian
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
#' Zurell, D., Franklin, J., König, C., Bouchet, P.J., Dormann, C.F., Elith, J.,
#' Fandos, G., Feng, X., Guillera‐Arroita, G., Guisan, A., Lahoz‐Monfort, J.J.,
#' Leitão, P.J., Park, D.S., Peterson, A.T., Rapacciuolo, G., Schmatz, D.R.,
#' Schröder, B., Serra‐Diaz, J.M., Thuiller, W., Yates, K.L., Zimmermann, N.E.
#' and Merow, C. (2020), A standard protocol for reporting species distribution
#' models. \emph{Ecography}, \strong{43}: 1261-1277.
#' @examples
#' \donttest{
#' library(raster)
#' # get predictor variables
#' f <- list.files(path=paste(system.file(package="phyloregion"),
#'                 '/ex', sep=''), pattern='.tif', full.names=TRUE)
#'
#' preds <- stack(f)
#' #plot(preds)
#' # get species occurrences
#' d <- read.csv(system.file("ex/Bombax.csv", package="phyloregion"))
#'
#' # fit ensemble model for four algorithms
#' mod <- sdm(d, predictors = preds)
#' }
#' @export
sdm <- function(x, pol = NULL, predictors = NULL, blank = NULL, tc = 2,
                lr = 0.001, bf = 0.75, n.trees = 50, step.size = n.trees, k=5,
                herbarium.rm = FALSE, n.points = 80, prob_method = "presab") {
    x <- .matchnames(x)
    name.sp <- unique(x$species)

    if (is.null(predictors)) {
        stop("you need to specify RasterStack of environmental data")
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
        herb_pol <- suppressWarnings(invisible(gBuffer(dx, width=0.5)))
        x <- x[is.na(over(x, herb_pol)),]
    }

    fam_pol <- pol

    if(nrow(x) < n.points){
        x <- .more_points(pts = x, preds = predictors)
        coordinates(x) <- ~lon+lat
        proj4string(x) <- CRS("+proj=longlat +datum=WGS84")
    }


    if (!is.null(pol)) {
        x <- x[pol,]
    }

    x1 <- x

    if (is.null(blank)) {
        blank <- predictors[[1]]
        blank[!is.na(blank)] <- 0
    }


    if (length(x1) > 0) {

        # OCCURRENCE POINTS
        occ_csv <- as.data.frame(x1)
        occ <- occ_csv[, c("lon", "lat")]

        set.seed(20220410)
        backg <- dismo::randomPoints(mask = predictors, n = nrow(occ)*3,
                                     ext = extent(x),
                                     extf = 1.25, warn = 0, p = occ)
        colnames(backg) <- c('lon', 'lat')

        # Arbitrarily assign group 1 as the testing data group
        test <- 1
        # Create vector of group memberships
        fold <- dismo::kfold(x = occ, k = k) # kfold is in dismo package

        # Separate observations into training and testing groups
        pres_train <- occ[fold != test, ]
        pres_test <- occ[fold == test, ]

        # Repeat the process for pseudo-absence points
        group <- dismo::kfold(x = backg, k = k)
        backg_train <- backg[group != test, ]
        backg_test <- backg[group == test, ]

        if (!is.null(fam_pol)) {
            predictors <- raster::crop(predictors, fam_pol)
        } else if (!is.null(pol)) {
            predictors <- raster::crop(predictors, pol)
        } else {
            predictors <- predictors
        }

        # only run if the maxent.jar file is available, in the right folder
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

        train <- rbind(pres_train, backg_train)
        pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
        envtrain <- extract(predictors, train)
        envtrain <- na.omit(data.frame(cbind(pa = pb_train, envtrain)))

        testpres <- data.frame(extract(predictors, pres_test))
        testbackg <- data.frame(extract(predictors, backg_test))

        Preds <- names(envtrain)[-1]
        y <- names(envtrain)[1]
        #Formula <- formula(paste(y," ~ ", paste(Preds, collapse=" + "),'+',
                           #paste("I(", Preds, "^2)", sep = "",
                                 #collapse = " + ")))

        Formula <- formula(paste(y," ~ ", paste(Preds, collapse=" + ")))

        if (prob_method == "presab") {
            # 1. Random Forest
        rf1 <- suppressWarnings(invisible(randomForest::randomForest(Formula,
                                                            data=envtrain)))
            erf <- dismo::evaluate(testpres, testbackg, rf1)
            px <- raster::predict(predictors, rf1,
                            ext = if (!is.null(pol)) extent(pol) else NULL)
            tr <- threshold(erf, 'spec_sens')
            trf <- px > tr
            # use the mask function
            zrf <- resample(trf, blank, method = "ngb")
            RF <- merge(zrf, blank)
            # TSS
            cm_rf <- data.frame(erf@confusion)
            tpr_rf <- cm_rf$tp / (cm_rf$tp + cm_rf$fn)
            tnr_rf <- cm_rf$tn / (cm_rf$fp + cm_rf$tn)
            tss_rf <- median(tpr_rf + tnr_rf - 1)
            gc()

            # 2. GLM
            gm <- suppressWarnings(invisible(glm(Formula,
                                        family = gaussian(link="identity"),
                                        data=envtrain)))
            ge <- dismo::evaluate(testpres, testbackg, gm)
            pg <- predict(predictors, gm,
                          ext = if (!is.null(pol)) extent(pol) else NULL)
            gtr <- threshold(ge, 'spec_sens')
            gt <- pg > gtr
            # use the mask function
            gt1 <- resample(gt, blank, method = "ngb")
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
                mx <- predict(predictors, xm,
                              ext = if (!is.null(pol)) extent(pol) else NULL,
                              progress='')
                mr <- threshold(xe, 'spec_sens')
                mt <- mx > mr
                # use the mask function
                mt1 <- resample(mt, blank, method = "ngb")
                MX <- merge(mt1, blank)
                # TSS
                cm_mx <- data.frame(xe@confusion)
                tpr_mx <- cm_mx$tp / (cm_mx$tp + cm_mx$fn)
                tnr_mx <- cm_mx$tn / (cm_mx$fp + cm_mx$tn)
                tss_mx <- median(tpr_mx + tnr_mx - 1)
            }
            gc()

            # 4. Generalised boosted models, GBMs
        gbm_mod <- suppressWarnings(invisible(dismo::gbm.step(data = envtrain,
                            gbm.x = Preds, gbm.y=y, family = "bernoulli",
                            tree.complexity = tc, learning.rate = lr,
                            bag.fraction = bf, verbose = TRUE, silent = TRUE,
                            plot.main = FALSE, step.size = step.size,
                            n.trees = n.trees)))

            egb <- suppressMessages(invisible(evaluate(testpres, testbackg,
                                                       gbm_mod)))
            pred_gb <- raster::predict(object = predictors, model = gbm_mod,
                                       n.trees = gbm_mod$gbm.call$best.trees,
                                       type = "response", na.rm = TRUE)
            gb <- gbm_mod$cv.statistics$cv.threshold
            gbf <- pred_gb > gb
            # use the mask function
            gbm1 <- resample(gbf, blank, method = "ngb")
            GBM <- merge(gbm1, blank)
            # TSS
            cm_gb <- data.frame(egb@confusion)
            tpr_gb <- cm_gb$tp / (cm_gb$tp + cm_gb$fn)
            tnr_gb <- cm_gb$tn / (cm_gb$fp + cm_gb$tn)
            tss_gbm <- median(tpr_gb + tnr_gb - 1)
            gc()

            if (maxent_available) {

                models <- list(RF, GLM, MX, GBM)
                models <- do.call(stack, models)
                names(models) <- c("RF", "GLM", "MAXENT", "GBM")
                m <- models[[which(!maxValue(models)==0)]]
                if (nlayers(m) > 1) {
                    m <- calc(m, median, forceapply=TRUE)
                }
                spo <- m==1

                aucs <- c(auc_RF=erf@auc, auc_GLM=ge@auc, auc_MAXENT=xe@auc,
                          auc_GBM=egb@auc)
                m_auc <- median(aucs)

                indiv_TSS <- c(TSS_RF=tss_rf, TSS_GLM=tss_glm,
                               TSS_MAXENT=tss_mx, TSS_GBM=tss_gbm)
                TSS <- median(indiv_TSS)

                occ_csv$species <- name.sp
                occ_csv$AUC <- m_auc
                occ_csv$TSS <- TSS
                occ_csv <- occ_csv[, c("species", "lon", "lat", "source",
                                       "AUC", "TSS")]
                occ_csv <- cbind(occ_csv, auc_RF=erf@auc, auc_GLM=ge@auc,
                                 auc_MAXENT=xe@auc, auc_GBM=egb@auc,
                                 TSS_RF=tss_rf, TSS_GLM=tss_glm,
                                 TSS_MAXENT=tss_mx, TSS_GBM=tss_gbm)
            } else {
                models <- list(RF, GLM, GBM)
                models <- do.call(stack, models)
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
                occ_csv <- occ_csv[, c("species", "lon", "lat", "source",
                                       "AUC", "TSS")]
                occ_csv <- cbind(occ_csv, auc_RF=erf@auc, auc_GLM=ge@auc,
                                 auc_GBM=egb@auc, TSS_RF=tss_rf,
                                 TSS_GLM=tss_glm, TSS_GBM=tss_gbm)
            }
        } else if (prob_method == "raw") {
            # 1. Random Forest
        rf1 <- suppressWarnings(invisible(randomForest::randomForest(Formula,
                                                            data=envtrain)))
            erf <- dismo::evaluate(testpres, testbackg, rf1)
            px <- raster::predict(predictors, rf1,
                            ext = if (!is.null(pol)) extent(pol) else NULL)
            #tr <- threshold(erf, 'spec_sens')
            #trf <- raster::crop(px, pol)
            # use the mask function
            zrf <- resample(px, blank, method = "ngb")
            #zrf <- resample(mask(trf, pol), blank, method = "ngb")
            RF <- merge(zrf, blank)
            # TSS
            cm_rf <- data.frame(erf@confusion)
            tpr_rf <- cm_rf$tp / (cm_rf$tp + cm_rf$fn)
            tnr_rf <- cm_rf$tn / (cm_rf$fp + cm_rf$tn)
            tss_rf <- median(tpr_rf + tnr_rf - 1)
            gc()

            # 2. GLM
            gm <- suppressWarnings(invisible(glm(Formula,
                                    family = gaussian(link="identity"),
                                    data=envtrain)))
            ge <- dismo::evaluate(testpres, testbackg, gm)
            pg <- predict(predictors, gm,
                          ext = if (!is.null(pol)) extent(pol) else NULL)
            #gtr <- threshold(ge, 'spec_sens')
            #gt <- raster::crop(pg, pol)
            # use the mask function
            gt1 <- resample(pg, blank, method = "ngb")
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
                mx <- predict(predictors, xm,
                              ext = if (!is.null(pol)) extent(pol) else NULL,
                              progress='')
                #mr <- threshold(xe, 'spec_sens')
                #mt <- raster::crop(mx, pol)
                # use the mask function
                mt1 <- resample(mx, blank, method = "ngb")
                MX <- merge(mt1, blank)
                # TSS
                cm_mx <- data.frame(xe@confusion)
                tpr_mx <- cm_mx$tp / (cm_mx$tp + cm_mx$fn)
                tnr_mx <- cm_mx$tn / (cm_mx$fp + cm_mx$tn)
                tss_mx <- median(tpr_mx + tnr_mx - 1)
            }
            gc()

            # 4. Generalised boosted models, GBMs
            gbm_mod <- suppressWarnings(invisible(dismo::gbm.step(data=envtrain,
                            gbm.x = Preds, gbm.y=y, family = "bernoulli",
                            tree.complexity = tc, learning.rate = lr,
                            bag.fraction = bf, verbose = TRUE, silent = TRUE,
                            plot.main = FALSE, step.size = step.size,
                            n.trees = n.trees)))

            egb <- suppressMessages(invisible(evaluate(testpres, testbackg,
                                                       gbm_mod)))
            pred_gb <- raster::predict(object = predictors, model = gbm_mod,
                                       n.trees = gbm_mod$gbm.call$best.trees,
                                       type = "response", na.rm = TRUE)
            #gb <- gbm_mod$cv.statistics$cv.threshold
            #gbf <- raster::crop(pred_gb, pol)
            # use the mask function
            gbm1 <- resample(pred_gb, blank, method = "ngb")
            GBM <- merge(gbm1, blank)
            # TSS
            cm_gb <- data.frame(egb@confusion)
            tpr_gb <- cm_gb$tp / (cm_gb$tp + cm_gb$fn)
            tnr_gb <- cm_gb$tn / (cm_gb$fp + cm_gb$tn)
            tss_gbm <- median(tpr_gb + tnr_gb - 1)
            gc()

            if (maxent_available) {
                models <- list(RF, GLM, MX, GBM)
                models <- do.call(stack, models)
                names(models) <- c("RF", "GLM", "MAXENT", "GBM")
                m <- models[[which(!maxValue(models)==0)]]
                if (nlayers(m) > 1) {
                    m <- calc(m, median, forceapply=TRUE)
                }
                spo <- m

                aucs <- c(auc_RF=erf@auc, auc_GLM=ge@auc, auc_MAXENT=xe@auc,
                          auc_GBM=egb@auc)
                m_auc <- median(aucs)

                indiv_TSS <- c(TSS_RF=tss_rf, TSS_GLM=tss_glm,
                               TSS_MAXENT=tss_mx, TSS_GBM=tss_gbm)
                TSS <- median(indiv_TSS)

                occ_csv$species <- name.sp
                occ_csv$AUC <- m_auc
                occ_csv$TSS <- TSS
                occ_csv <- occ_csv[, c("species", "lon", "lat", "source",
                                       "AUC", "TSS")]
                occ_csv <- cbind(occ_csv, auc_RF=erf@auc, auc_GLM=ge@auc,
                                 auc_MAXENT=xe@auc, auc_GBM=egb@auc,
                                 TSS_RF=tss_rf, TSS_GLM=tss_glm,
                                 TSS_MAXENT=tss_mx, TSS_GBM=tss_gbm)
            } else {
                models <- list(RF, GLM, GBM)
                models <- do.call(stack, models)
                names(models) <- c("RF", "GLM", "GBM")

                m <- models[[which(!maxValue(models)==0)]]
                if (nlayers(m) > 1) {
                    m <- calc(m, median, forceapply=TRUE)
                }
                spo <- m

                aucs <- c(auc_RF=erf@auc, auc_GLM=ge@auc, auc_GBM=egb@auc)
                m_auc <- median(aucs)

                indiv_TSS <- c(TSS_RF=tss_rf, TSS_GLM=tss_glm, TSS_GBM=tss_gbm)
                TSS <- median(indiv_TSS)

                occ_csv$species <- name.sp
                occ_csv$AUC <- m_auc
                occ_csv$TSS <- TSS
                occ_csv <- occ_csv[, c("species", "lon", "lat", "source",
                                       "AUC", "TSS")]
                occ_csv <- cbind(occ_csv, auc_RF=erf@auc, auc_GLM=ge@auc,
                                 auc_GBM=egb@auc, TSS_RF=tss_rf,
                                 TSS_GLM=tss_glm, TSS_GBM=tss_gbm)
            }
        }

        return(list(ensemble_raster = spo, ensemble_AUC = m_auc,
                    ensemble_TSS=TSS, data=occ_csv,
                    indiv_models=models, indiv_AUCs=aucs, indiv_TSS=indiv_TSS))
        gc()
    }
}





