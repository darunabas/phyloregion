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

myRF <- function(x, p, q, blank) {
    tf <- tuneRF(x[, 2:ncol(x)], x[, "pa"], plot=FALSE)
    mt <- tf[which.min(tf[,2]), 1]
    rf1 <- randomForest(x[, 2:ncol(x)], x[, "pa"], mtry=mt, ntree=250,
                        na.action = na.omit)
    rp <- predict(p, rf1, na.rm=TRUE)
    erf <- pa_evaluate(predict(rf1, q[q$pa==1, ]), predict(rf1, q[q$pa==0, ]))
    tr <- erf@thresholds
    trf <- rp > tr$max_spec_sens
    zrf <- terra::resample(trf, blank, method = "near")
    rf_fnl <- merge(zrf, blank)
    return(RF=rf_fnl)
}

myGLM <- function(Formula, x, p, q, blank) {
    gm <- glm(Formula, family = gaussian(link = "identity"), data = x)
    pg <- predict(p, gm, type="response")
    ge <- pa_evaluate(predict(gm, q[q$pa==1, ]), predict(gm, q[q$pa==0, ]))
    gtr <- ge@thresholds
    gt <- pg > gtr$max_spec_sens
    gt1 <- terra::resample(gt, blank, method = "near")
    glm_fnl <- merge(gt1, blank)
    return(GLM=glm_fnl)
}

myMAXENT <- function(p, ox, ot, blank, xy) {
    maxent_available <- FALSE
    if (maxentropy()) {
        maxent_available <- TRUE
        xm <- maxentropy(p, ox)
        mx <- predict(xm, p)
        xe <- pa_evaluate(xm, p = ot, a = xy, x = p)
        mr <- xe@thresholds
        mt <- mx > mr$max_spec_sens
        mt1 <- terra::resample(mt, blank, method = "near")
        MX <- merge(mt1, blank)
        return(MAXENT=MX)
    }
}

# write functions for each model
# model=c("RF", "GLM", "MAXENT", "GBM") to choose models
# run bioclim to generate points for species with limited sampling.

sampleBuffer <- function(p, size, width) {
    if (max(distance(p)) > 7000000) {
        h <- buffer(p, width = width)
    } else {
        h <- convHull(p)
    }
    spatSample(h, size = size, method="random")
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
#' @param algorithm Character. The choice of algorithm to run the species
#' distribution model. Available algorithms include:
#' \itemize{
#' \item \dQuote{RF}: Random forest.
#' \item \dQuote{GLM}: Generalized linear model.
#' \item \dQuote{MAXENT}: Maximum entropy.
#' \item \dQuote{GBM}: Generalized boosted regressions model.}
#'
#' @param size Minimum number of points required to successfully run
#' a species distribution model
#' @param width Width of buffer in meter if x is in longitude/latitude CRS.
#' @rdname sdmterra
#' @importFrom terra distance convHull spatSample vect ext window<- rast nlyr
#' @importFrom terra geom resample
#' @importFrom predicts make_folds maxentropy pa_evaluate
#' @importFrom randomForest randomForest tuneRF
#' @importFrom stats glm median formula gaussian
#' @importFrom grDevices chull
#' @importFrom gbm gbm
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
#' @export
sdmterra <- function(x, predictors = NULL, blank = NULL,
                     algorithm = "GLM",
                     size = 50, width = 100000) {
    x <- .matchnames(x)
    name.sp <- unique(x$species)

    if (is.null(predictors)) {
        stop("you need to specify RasterStack of environmental data")
    }

    if (is.null(blank)) {
        blank <- predictors[[1]]
        blank[!is.na(blank)] <- 0
    }

    x <- x[, -1]
    x <- unique(na.omit(x))
    x$source <- "raw"
    x <- vect(x, crs="+proj=longlat +datum=WGS84")

    if (length(x) < size) {
        h <- geom(x, df = TRUE)
        h <- h[, 3:4]
        h$source <- "raw"
        v <- sampleBuffer(x, size = size-length(x), width = width)
        v <- geom(v, df = TRUE)
        v <- v[, 3:4]
        v$source <- "random"
        x <- rbind(h, v)
    } else {
        x <- geom(x, df = TRUE)
        x <- x[, 3:4]
        x$source <- "raw"
    }

    occ_csv <- x
    x <- x[, -3]

    b <- extract(predictors, x[, 1:2])
    b <- b[, -1]

    #set.seed(0)
    e <- ext(vect(x, geom = c("x", "y"))) + 1
    window(predictors) <- e
    bg <- spatSample(predictors, size = nrow(x)*100, "random",
                     na.rm = TRUE, xy = TRUE)
    xy <- bg[, 1:2]
    bg <- bg[, -c(1:2)]
    j <- data.frame(rbind(cbind(pa = 1, b), cbind(pa = 0, bg)))
    j <- na.omit(j)

    # for random forest
    i <- sample(nrow(j), 0.2 * nrow(j))
    test <- j[i,]
    train <- j[-i,]

    # for maxent
    fold <- predicts::make_folds(x, k=5)
    occtest <- x[fold == 1, ]
    occtrain <- x[fold != 1, ]

    # for glm
    Preds <- names(train)[-1]
    y <- names(train)[1]
    Formula <- formula(paste(y," ~ ", paste(Preds, collapse=" + ")))

    if (algorithm == "all") {
        models <- list()
        # 1. Random Forest
        tryCatch({
            RF1 <- FALSE
            rf_fnl <- myRF(x = train, p = predictors, q=test,
                 blank=blank)
            RF1 <- TRUE
            if (RF1) models[[1]] <- rf_fnl
            names(models[[1]]) <- "RF"
        }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        # 2. GLM
        tryCatch({
            GLM1 <- FALSE
            glm_fnl <- myGLM(Formula = Formula, x = train, p = predictors,
                             q = test, blank = blank)
            GLM1 <- TRUE
            if (GLM1) models[[2]] <- glm_fnl
            names(models[[2]]) <- "GLM"
        }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        # 3. Maxent
        tryCatch({
            MX1 <- FALSE
            MX <- myMAXENT(p = predictors, ox = occtrain, xy = xy,
                           ot = occtest, blank = blank)
            MX1 <- TRUE
            if (MX1) models[[3]] <- MX
            names(models[[3]]) <- "MX"
        }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        models <- rast(models)

    } else {
        models <- switch(algorithm,
                         RF = myRF(x = train, p = predictors, q=test,
                                   blank=blank),
                         GLM = myGLM(Formula = Formula, x = train, p = predictors,
                                     q = test, blank = blank),
                         MAXENT = myMAXENT(p = predictors, ox = occtrain, xy = xy,
                                           ot = occtest, blank = blank))
        names(models) <- algorithm
    }

    if (nlyr(models) > 1) {
        models <- median(models)
    }
    spo <- models > 0.5
    occ_csv$species <- name.sp

    return(list(ensemble_raster = spo, data=occ_csv, indiv_models=models))

}





