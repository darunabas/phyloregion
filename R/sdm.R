boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                  rm.duplicate = TRUE, method = 'spearman') {
  
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    return(round(pi/ei,10))
  }
  
  if (inherits(fit,"SpatRaster")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- extract(fit, obs, ID=FALSE, raw=TRUE) |> na.omit()
    }
    fit <- values(fit)
    fit <- fit[!is.na(fit)]
  }
  
  mini <- min(fit,obs)
  maxi <- max(fit,obs)
  
  if(length(nclass)==1){
    if (nclass == 0) { #moving window
      if (window.w == "default") {window.w <- (max(fit) - min(fit))/10}
      vec.mov <- seq(from = mini, to = maxi - window.w, 
                     by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else{ #window based on nb of class
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else{ #user defined window
    vec.mov <- c(mini, sort(nclass[!nclass>maxi|nclass<mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)  
  }
  return(boyce = round(b, 3))
}


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

myGLM <- function(f, x, p, q, forcast) {
  gm <- glm(f, family = gaussian(link = "identity"), data = x)
  pg <- terra::predict(p, gm, type="response")
  ge <- predicts::pa_evaluate(predict(gm, q[q$pa==1, ]), 
                              predict(gm, q[q$pa==0, ]))
  gtr <- ge@thresholds
  
  pg <- terra::predict(forcast, gm, type="response")
  glm_fnl <- pg > gtr$max_spec_sens
  
  # TSS
  cm_gl <- data.frame(ge@confusion)
  tpr_gl <- cm_gl$tp / (cm_gl$tp + cm_gl$fn)
  tnr_gl <- cm_gl$tn / (cm_gl$fp + cm_gl$tn)
  tss_glm <- median(tpr_gl + tnr_gl - 1)
  
  # AUC
  auc_glm <- ge@stats$auc
  df_glm <- data.frame(TSS = tss_glm, AUC = auc_glm, thresh=gtr$max_spec_sens)
  return(list(GLM = glm_fnl, data_glm = df_glm, GLM_raw = pg))
}

myMAXENT <- function(p, ox, ot, xy, forcast) {
  maxent_available <- FALSE
  if (MaxEnt()) {
    maxent_available <- TRUE
    performanceMaxent <- numeric(0)
    for( beta in c(2, 5, 10, 15, 20, 30)) {
      xm <- MaxEnt(p, ox, args = c(
        paste0('betamultiplier=', beta),
        'linear=true',
        'quadratic=false',
        'product=false',
        'threshold=true',
        'hinge=true',
        'responsecurves=false',
        'jackknife=false'
      ))
      mx <- terra::predict(xm, p)
      xe <- predicts::pa_evaluate(xm, p = ot, a = xy, x = p)
      performanceMaxent <- c(performanceMaxent, xe@stats$auc)
    }
    
    beta <- c(2, 5, 10, 15, 20, 30)[which.max(performanceMaxent)]
    xm <- MaxEnt(p, ox, args = c(
      paste0('betamultiplier=', beta),
      'linear=true',
      'quadratic=false',
      'product=false',
      'threshold=true',
      'hinge=true',
      'responsecurves=false',
      'jackknife=false'
    ))
    
    mx <- terra::predict(xm, p)
    xe <- predicts::pa_evaluate(xm, p = ot, a = xy, x = p)
    mr <- xe@thresholds
    
    mx <- terra::predict(xm, forcast)
    MX <- mx > mr$max_spec_sens
    
    # TSS
    cm_mx <- data.frame(xe@confusion)
    tpr_mx <- cm_mx$tp / (cm_mx$tp + cm_mx$fn)
    tnr_mx <- cm_mx$tn / (cm_mx$tn + cm_mx$fp)
    tss <- median(tpr_mx + tnr_mx - 1)
    
    # AUC
    aucs <- xe@stats$auc
    stats_df <- data.frame(TSS = tss, AUC = aucs, thresh=mr$max_spec_sens)
    return(list(MAXENT = MX, data_maxent = stats_df, MAXENT_raw = mx))
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

#' Fast species distribution model
#'
#' This function computes species distribution models using
#' two modelling algorithms: generalized linear models,
#' and maximum entropy (only if \code{rJava} is available).
#' Note: this is an experimental function, and may change in the future.
#'
#' @param x A dataframe containing the species occurrences
#' and geographic coordinates. Column 1 labeled as "species", column 2 "lon",
#' column 3 "lat".
#' @param layers A \code{SpatRaster} of predictor variables for fitting species
#' distribution models from species occurrences.
#' @param predictors If predicting to new time points, the climate layers for
#' the time points.
#' @param background A dataframe of background points, specifying 2 columns 
#' with long lat or x and y as nulls for species distribution modeling, often 
#' using a vector of probability weights.
#' @param pol A vector polygon specifying the calibration area or boundary to 
#' account for a more realistic dispersal capacity and ecological limitation 
#' of a species. If \code{NULL}, the extent of input points is used.
#' @param algorithm Character. The choice of algorithm to run the species
#' distribution model. For now, the available algorithms include:
#' \itemize{
#' \item \dQuote{all}: Calls all available algorithms: both GLM and MAXENT.
#' \item \dQuote{GLM}: Calls only Generalized linear model.
#' \item \dQuote{MAXENT}: Calls only Maximum entropy.
#' }
#' @param size Minimum number of points required to successfully run
#' a species distribution model especially for species with few occurrences.
#' @param thin Whether to spatially thin occurrences
#' @param thin.size The size of the thin occurrences.
#' @param width Width of buffer in meter if x is in longitude/latitude CRS.
#' @param mask logical. Should \code{layers} be used to mask? Only used if 
#' \code{pol} is a SpatVector.
#' @rdname sdm
#' @importFrom terra distance convHull spatSample vect ext window<- rast nlyr
#' @importFrom terra geom resample crop median deepcopy as.polygons predict
#' @importFrom terra extract values
#' @importFrom predicts folds MaxEnt pa_evaluate
#' @importFrom stats glm median formula gaussian dist runif reformulate cor
#' @importFrom smoothr smooth
#' @return A list with the following objects:
#' \itemize{
#'   \item \code{ensemble_raster} The ensembled raster that predicts
#'   the potential species distribution based on the algorithms selected.
#'   \item \code{data} Model evaluation statistics: true skill statistic
#'   (\code{TSS}), area under the curve (\code{AUC}) and the threshold that
#'   maximizes sensitivity plus specificity (\code{thresh}); for a single
#'   algorithm also the Boyce index (\code{boyce}). With
#'   \code{algorithm = "all"}, the median across algorithms.
#'   \item \code{polygon} Map polygons of the predicted distributions
#'   analogous to extent-of-occurrence range polygon.
#'   \item \code{indiv_models} Raster layers for the separate models that
#'   predict the potential species distribution.
#'   \item \code{raw_habitatsuit} Continuous habitat suitability before
#'   thresholding, one layer per algorithm.
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
#' # m <- sdm(d, layers = preds, predictors = preds, algorithm = "all")
#' # plot(m$ensemble_raster)
#' # plot(m$polygon, add=TRUE)
#' }
#' @export
sdm <- function (x, layers = NULL, pol = NULL, thin = TRUE, thin.size = 500,
                   algorithm = "all", size = 50, width = 50000, 
                   mask = FALSE, predictors, background = NULL) {
  
  x <- .matchnames(x)
  name.sp <- unique(x$species)
  
  if (is.null(layers)) {
    stop("you need to specify raster layers of environmental data")
  }
  
  x <- x[, -1]
  x <- unique(na.omit(x))
  
  if (thin == TRUE) {
    x <- thin_max(x = x, cols = c("lat", "lon"), npoints = thin.size)
  }
  x$source <- "rw"
  x <- vect(x, crs = "EPSG:4326")
  if (length(x) < size) {
    h <- geom(x, df = TRUE)
    h <- h[, 3:4]
    h$source <- "rw"
    v <- sampleBuffer(x, size = size - length(x), width = width)
    v <- geom(v, df = TRUE)
    v <- v[, 3:4]
    v$source <- "rn"
    x <- rbind(h, v)
  }
  else {
    x <- geom(x, df = TRUE)
    x <- x[, 3:4]
    x$source <- "rw"
  }
  
  x <- x[, -3]
  b <- extract(layers, x[, 1:2], ID = FALSE)
  
  if(is.null(pol)) {
    pol <- ext(vect(x, geom = c("x", "y"))) + 1
  }
  else {
    pol <- pol
  }

  layers <- terra::crop(layers, pol, mask = mask)
  predictors <- terra::crop(predictors, pol, mask = mask)
  
  if(is.null(background)) {
    background <- spatSample(layers, size = nrow(x)*100, method = "random",
                     na.rm = TRUE, xy = TRUE, warn = FALSE, ext=ext(pol))
    xy <- background[, 1:2]
  } else {
    xy <- background[, 1:2]
  }
  
  bg <- terra::extract(layers, xy, ID = FALSE)
  
  j <- data.frame(rbind(cbind(pa = 1, b), cbind(pa = 0, bg))) |> na.omit()
  i <- sample(nrow(j), 0.25 * nrow(j))
  
  test <- j[i, ]
  train <- j[-i, ]
  fold <- predicts::folds(x, k = 5)
  occtest <- x[fold == 1, ]
  occtrain <- x[fold != 1, ]

  Formula <- reformulate(termlabels = colnames(train)[-1], response="pa")

  if (algorithm == "all") {
    models <- c()
    raws <- list()
    dat <- list()
    tryCatch({
      GLM1 <- FALSE
      glm_fnl <- myGLM(f = Formula, x = train, p = layers,
                       q = test, forcast = predictors)
      GLM1 <- TRUE
      if (GLM1) {
        names(glm_fnl)[1] <- "GLM"
        models <- c(models, glm_fnl[[1]])
        raws$GLM <- glm_fnl[[3]]
        dat[[1]] <- glm_fnl[[2]]
      }
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
    })
    tryCatch({
      MX1 <- FALSE
      MX <- myMAXENT(p = layers, ox = occtrain, xy = xy,
                     ot = occtest, forcast = predictors)
      MX1 <- TRUE
      if (MX1) {
        names(MX)[[1]] <- "MAXENT"
        models <- c(models, MX[[1]])
        raws$MAXENT <- MX[[3]]
        dat[[2]] <- MX[[2]]
      }
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
    })
    models <- rast(models)
    raw_model <- rast(raws)
    names(raw_model) <- names(raws)
    dat <- do.call(rbind, dat)
    dat <- apply(dat, 2, FUN=median)
  }
  else {
    models <- switch(algorithm,
                     GLM = myGLM(f = Formula, x = train, p = layers,
                                 q = test, forcast = predictors),
                     MAXENT = myMAXENT(p = layers, ox = occtrain, xy = xy,
                                       ot = occtest, forcast = predictors))
    names(models)[[1]] <- algorithm
    dat <- models[[2]]
    boy <- boyce(fit = models[[3]], obs = occtest)
    raw_model <- models[[3]]
    dat <- data.frame(dat, boyce = boy)
    models <- models[[1]]
    
  }
  m <- models
  if (nlyr(m) > 1) {
    m <- terra::median(m)
  }
  spo <- m > 0.5
  pol <- smoothr::smooth(as.polygons(spo), method = "ksmooth")
  names(pol)[1] <- "median"
  pol <- pol[pol$median == 1, ]
  pol$species <- name.sp
 
  return(list(ensemble_raster = spo, data = dat, polygon = pol,
              indiv_models = models, raw_habitatsuit = raw_model))
}




