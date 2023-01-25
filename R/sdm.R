thin_max <- function(x, cols, npoints){
    inds <- vector(mode="numeric")
    this.dist <- as.matrix(dist(x[,cols], upper=TRUE))
    inds <- c(inds, as.integer(runif(1, 1, length(this.dist[,1]))))
    inds <- c(inds, which.max(this.dist[,inds]))
    while(length(inds) < npoints){
        min.dists <- apply(this.dist[,inds], 1, min)
        this.ind <- which.max(min.dists)
        if(length(this.ind) > 1){
            print("Breaking tie...")
            this.ind <- this.ind[1]
        }
        inds <- c(inds, this.ind)
    }
    return(x[inds,])
}

#myRF <- function(x, p, q) {
#    tf <- tuneRF(x[, 2:ncol(x)], x[, "pa"], plot=FALSE)
#    mt <- tf[which.min(tf[,2]), 1]
#    rf1 <- randomForest(x[, 2:ncol(x)], x[, "pa"], mtry=mt, ntree=250,
#                        na.action = na.omit)
#    rp <- terra::predict(p, rf1, na.rm=TRUE)
#    erf <- pa_evaluate(predict(rf1, q[q$pa==1, ]), predict(rf1, q[q$pa==0, ]))
#    tr <- erf@thresholds
#    rf_fnl <- rp > tr$max_spec_sens
#    return(RF = rf_fnl)
#}

myGLM <- function(f, x, p, q) {
    gm <- glm(f, family = gaussian(link = "identity"), data = x)
    pg <- terra::predict(p, gm, type="response")
    ge <- pa_evaluate(predict(gm, q[q$pa==1, ]), predict(gm, q[q$pa==0, ]))
    gtr <- ge@thresholds
    glm_fnl <- pg > gtr$max_spec_sens
    return(GLM = glm_fnl)
}

myMAXENT <- function(p, ox, ot, xy) {
    maxent_available <- FALSE
    if (MaxEnt()) {
        maxent_available <- TRUE
        xm <- MaxEnt(p, ox)
        mx <- terra::predict(xm, p)
        xe <- pa_evaluate(xm, p = ot, a = xy, x = p)
        mr <- xe@thresholds
        MX <- mx > mr$max_spec_sens
        return(MAXENT=MX)
    }
}

#myCTA <- function(f, x, p, q) {
#    cm <- rpart(f, data = x)
#    pc <- terra::predict(p, cm)
#    ce <- pa_evaluate(predict(cm, q[q$pa==1, ]), predict(cm, q[q$pa==0, ]))
#    ctr <- ce@thresholds
#    cta_fnl <- pc > ctr$max_spec_sens
#    return(CTA = cta_fnl)
#}

#myGBM <- function(f, x, p, q) {
#    gbm_mod <- gbm::gbm(f, data = x, interaction.depth = 3,
#                        n.minobsinnode = 1, cv.folds = 5)
#
#    pred_gb <- terra::predict(p, gbm_mod, na.rm=TRUE)
#    egb <- predicts::pa_evaluate(predict(gbm_mod, q[q$pa==1, ]),
#                                 predict(gbm_mod, q[q$pa==0, ]))
#    gb <- egb@thresholds
#    gbm_fnl <- pred_gb > gb$max_spec_sens
#    return(GBM = gbm_fnl)
#}

sampleBuffer <- function(p, size, width) {
    if (max(distance(p)) > 7000000) {
        h <- buffer(p, width = width)
    } else {
        h <- convHull(p)
    }
    spatSample(h, size = size, method = "random")
}

#' Species distribution models
#'
#' This function computes species distribution models using
#' two modelling algorithms: generalized linear models,
#' and maximum entropy (only if \code{rJava} is available).
#' Note: this is an experimental function, and may change in the future.
#'
#' @param x A dataframe containing the species occurrences
#' and geographic coordinates. Column 1 labeled as "species", column 2 "lon",
#' column 3 "lat".
#' @param predictors A \code{SpatRaster} to extract values from the
#' locations in x on which the models will be projected.
#' @param pol A vector polygon specifying the boundary to restrict the
#' prediction. If \code{NULL}, the extent of input points is used.
#' @param algorithm Character. The choice of algorithm to run the species
#' distribution model. Available algorithms include:
#' \itemize{
#' \item \dQuote{all}: Calls all available algorithms: GLM, and MAXENT.
#' \item \dQuote{GLM}: Calls only Generalized linear model.
#' \item \dQuote{MAXENT}: Calls only Maximum entropy.
#' }
#' @param size Minimum number of points required to successfully run
#' a species distribution model especially for species with few occurrences.
#' @param thin Whether to thin occurrences
#' @param thin.size The size of the thin occurrences.
#' @param width Width of buffer in meter if x is in longitude/latitude CRS.
#' @param mask logical. Should y be used to mask? Only used if pol is a SpatVector
#' @rdname sdm
#' @importFrom terra distance convHull spatSample vect ext window<- rast nlyr
#' @importFrom terra geom resample crop median deepcopy as.polygons predict
#' @importFrom terra extract
#' @importFrom predicts folds MaxEnt pa_evaluate
#' @importFrom stats glm median formula gaussian dist runif
#' @importFrom smoothr smooth
#' @return A list with the following objects:
#' \itemize{
#'   \item \code{ensemble_raster} The ensembled raster that predicts
#'   the potential species distribution based on the algorithms selected.
#'   \item \code{data} The dataframe of occurrences used to implement the model.
#'   \item \code{polygon} Map polygons of the predicted distributions
#'   analogous to extent-of-occurrence range polygon.
#'   \item \code{indiv_models} Raster layers for the separate models that
#'   predict the potential species distribution.
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
#' # get predictor variables
#' library(predicts)
#' f <- system.file("ex/bio.tif", package="predicts")
#' preds <- rast(f)
#' #plot(preds)
#'
#' # get species occurrences
#' b <- file.path(system.file(package="predicts"), "ex/bradypus.csv")
#' d <- read.csv(b)
#'
#' # fit ensemble model for four algorithms
#' m <- sdm(d, predictors = preds, algorithm = "all")
#' # plot(m$ensemble_raster)
#' # plot(m$polygon, add=TRUE)
#' }
#' @export
sdm <- function(x, predictors = NULL, pol = NULL, thin = TRUE, thin.size = 500,
                algorithm = "all", size = 50, width = 50000, mask = FALSE) {
    x <- .matchnames(x)
    name.sp <- unique(x$species)

    if (is.null(predictors)) {
        stop("you need to specify RasterStack of environmental data")
    }

    x <- x[, -1]
    x <- unique(na.omit(x))

    if(thin == TRUE) {
        x <- thin_max(x=x, cols = c("lat", "lon"), npoints = thin.size)
    }
    x$source <- "rw"
    x <- vect(x, crs="+proj=longlat +datum=WGS84")

    if (length(x) < size) {
        h <- geom(x, df = TRUE)
        h <- h[, 3:4]
        h$source <- "rw"
        v <- sampleBuffer(x, size = size-length(x), width = width)
        v <- geom(v, df = TRUE)
        v <- v[, 3:4]
        v$source <- "rn"
        x <- rbind(h, v)
    } else {
        x <- geom(x, df = TRUE)
        x <- x[, 3:4]
        x$source <- "rw"
    }

    occ_csv <- x
    x <- x[, -3]

    b <- extract(predictors, x[, 1:2])
    b <- b[, -1]

    if(!is.null(pol)) {
        pol <- buffer(pol, width=width)
    } else {
        pol <- ext(vect(x, geom = c("x", "y"))) + 1
    }

    predictors <- terra::crop(predictors, pol, mask = mask)
    bg <- spatSample(predictors, size = nrow(x)*100, method = "random",
                     na.rm = TRUE, xy = TRUE, warn = FALSE, ext=ext(pol))
    xy <- bg[, 1:2]
    bg <- bg[, -c(1:2)]
    j <- data.frame(rbind(cbind(pa = 1, b), cbind(pa = 0, bg)))
    j <- na.omit(j)

    # for random forest
    i <- sample(nrow(j), 0.25 * nrow(j))
    test <- j[i,]
    train <- j[-i,]

    # for maxent
    fold <- predicts::folds(x, k=5)
    occtest <- x[fold == 1, ]
    occtrain <- x[fold != 1, ]

    # for glm
    Preds <- names(train)[-1]
    y <- names(train)[1]
    Formula <- formula(paste(y," ~ ", paste(Preds, collapse=" + ")))

    if (algorithm == "all") {
        models <- c()

        # 1. Random Forest
        #tryCatch({
        #    RF1 <- FALSE
        #    rf_fnl <- myRF(x = train, p = predictors, q = test)
        #    RF1 <- TRUE
        #    if (RF1) {
        #        names(rf_fnl) <- "RF"
        #        models <- c(models, rf_fnl)
        #    }
        #}, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        # 2. GLM
        tryCatch({
            GLM1 <- FALSE
            glm_fnl <- myGLM(f = Formula, x = train, p = predictors, q = test)
            GLM1 <- TRUE
            if (GLM1) {
                names(glm_fnl) <- "GLM"
                models <- c(models, glm_fnl)
            }
        }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        # 3. Maxent
        tryCatch({
            MX1 <- FALSE
            MX <- myMAXENT(p=predictors, ox=occtrain, xy=xy, ot=occtest)
            MX1 <- TRUE
            if (MX1) {
                names(MX) <- "MAXENT"
                models <- c(models, MX)
            }
        }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        # 4. gbm
        #tryCatch({
        #    GBM1 <- FALSE
        #    gbm_fnl <- myGBM(f = Formula, x = train, p = predictors, q = test)
        #    GBM1 <- TRUE
        #    if (GBM1) {
        #        names(gbm_fnl) <- "GBM"
        #        models <- c(models, gbm_fnl)
        #    }
        #}, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        # 5. CTA
        #tryCatch({
        #    CTA1 <- FALSE
        #    cta_fnl <- myCTA(f = Formula, x = train, p = predictors, q = test)
        #    CTA1 <- TRUE
        #    if (CTA1) {
        #        names(cta_fnl) <- "CTA"
        #        models <- c(models, cta_fnl)
        #    }
        #}, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})

        models <- rast(models)
    } else {
        models <- switch(algorithm,
                         #RF = myRF(x = train, p = predictors, q = test),
                         GLM = myGLM(f=Formula, x=train, p=predictors, q=test),
                         MAXENT = myMAXENT(p=predictors, ox=occtrain, xy=xy,
                                           ot=occtest),
                         #GBM = myGBM(f=Formula, x=train, p=predictors, q=test),
                         #CTA = myCTA(f=Formula, x=train, p=predictors, q=test)
                         )
        names(models) <- algorithm
    }
    m <- models
    if (nlyr(m) > 1) {
        m <- terra::median(m)
    }
    spo <- m > 0.5
    pol <- smooth(as.polygons(spo), method = "ksmooth")
    names(pol)[1] <- "median"
    pol <- pol[pol$median==1,]
    pol$species <- name.sp
    occ_csv$species <- name.sp
    return(list(ensemble_raster = spo, data = occ_csv,
                polygon = pol, indiv_models = models))
}





