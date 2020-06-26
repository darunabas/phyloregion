matching <- function (x, y)
{
  if (ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  x <- data.matrix(x)
  y <- data.matrix(y)
  z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (i in 1:nrow(y)) {
    z[, i] <- colSums(1/((t(x) - y[i, ])/(t(x) - y[i, ])),
                      na.rm = TRUE)/ncol(x)
  }
  rownames(z) <- rownames(x)
  colnames(z) <- rownames(y)
  return(z)
}


cooccur <- function (data)
{
  if ((is.matrix(data) || is.data.frame(data)) == FALSE)
    stop("The data must be a matrix or a data frame!")
  rn <- rownames(data)
  if (is.logical(data) == TRUE)
    data <- data.matrix(as.data.frame(data)) + as.integer(1)
  if (is.numeric(data) == TRUE)
    data <- apply(data, 2, function(x) as.integer(x))
  if (is.data.frame(data) == TRUE)
    data <- apply(data, 2, function(x) as.integer(as.factor(x)))
  col <- ncol(data)
  if (col != 1) {
    newdist <- function(data, col = col, colnum) {
      nvar <- 1:col
      n <- length(levels(as.factor(data[, colnum])))
      var <- prob <- vector("list", (col - 1))
      for (i in 1:(col - 1)) {
        var[[i]] <- table(data[, nvar[-colnum][i]],
                          data[, colnum])
        prob[[i]] <- var[[i]]/matrix(colSums(var[[i]]),
                                     nrow = nrow(var[[i]]), ncol = ncol(var[[i]]),
                                     byrow = TRUE)
      }
      probmat <- do.call(rbind, prob)
      matnew <- matrix(0, nrow = n, ncol = n)
      rownames(matnew) <- colnames(matnew) <- 1:n
      for (i in 1:n) {
        for (j in 1:n) {
          matnew[i, j] <- (sum(apply(probmat[, c(i,
                                                 j)], 1, max)) - (col - 1))/(col - 1)
        }
      }
      return(matnew)
    }
    newdata <- vector("list", col)
    for (i in 1:col) {
      newdata[[i]] <- newdist(data, col = col, i)
    }
    distmat <- matrix(0, nrow(data), nrow(data))
    for (i in 1:nrow(data)) {
      for (j in 1:nrow(data)) {
        distsum <- numeric(col)
        for (k in 1:col) {
          distsum[k] <- newdata[[k]][data[i, k], data[j,
                                                      k]]
        }
        distmat[i, j] <- sum(distsum)
      }
    }
    rownames(distmat) <- colnames(distmat) <- rn
  }
  else {
    distmat <- matching(data, data)
    rownames(distmat) <- colnames(distmat) <- rn
    warning("Due to only 1 variable, simple matching distance is calculated instead!\n            To produce coocurrence distance, it requires at least 2 variables.")
  }
  return(distmat)
}




distNumeric <- function (x, y, method = "mrw", xyequal = TRUE)
{
  if ((is.matrix(x) && is.matrix(x)) == FALSE)
    stop("x and y must be a matrix object!")
  if (ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  if (xyequal == TRUE) {
    span <- apply(x, 2, function(x) max(x) - min(x))
    variance <- apply(x, 2, var)
  }
  else {
    span <- apply(rbind(x, y), 2, function(x) max(x) - min(x))
    variance <- apply(rbind(x, y), 2, var)
  }
  span.2 <- span^2
  num_distance <- c("mrw", "sev", "ser", "ser.2", "se")
  method <- match.arg(method, num_distance)
  result <- switch(method,
                   mrw = weightedNum(x, y, p = 1, alpha = span),
                   ser = weightedNum(x, y, p = 2, alpha = span),
                   ser.2 = weightedNum(x, y, p = 2, alpha = span.2),
                   sev = weightedNum(x, y, p = 2, alpha = variance),
                   se = weightedNum(x, y, p = 2, alpha = 1))
  return(result)
}


weightedNum <- function (x, y, p = 2, alpha = 1)
{
  if (ncol(x) != ncol(y))
    stop(sQuote("x"), " and ", sQuote("y"), " must have the same number of columns")
  x <- data.matrix(x)
  y <- data.matrix(y)
  z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (i in 1:nrow(y)) {
    z[, i] <- colSums(1/alpha * abs(t(x) - y[i, ])^p, na.rm = TRUE)
  }
  rownames(z) <- rownames(x)
  colnames(z) <- rownames(y)
  return(z)
}

gower_distance <- function (data, method = "gower", idnum = NULL, idbin = NULL,
                            idcat = NULL) {
  if (any(is.na(data)))
    stop("Cannot handle missing values!")
  if ((is.matrix(data) || is.data.frame(data)) == FALSE)
    stop("The data must be a matrix or a data frame object!")
  if (is.null(idnum) && is.null(idbin) && is.null(idcat))
    stop("There is no distance computation, specify the numerical, binary, categorical variables!")
  if (is.null(idbin) && is.null(idcat) || is.null(idnum) &&
      is.null(idcat) || is.null(idnum) && is.null(idbin))
    stop("There is no mixed variables!")
  dist_num4 <- c("gower", "wishart", "podani", "huang", "harikumar",
                 "ahmad")
  method <- match.arg(method, dist_num4)
  if ((length(c(idbin, idcat)) == 1) & method == "ahmad")
    stop("Ahmad-Dey distance can not be calculated\n         because the combined binary and categorical variable is only 1 variable!")
  if (length(idcat) == 1 & method == "harikumar")
    stop("Harikumar-PV distance can not be calculated\n         because the categorical variable is only 1 variable!")
  if (is.null(idnum)) {
    num <- 0
    msd <- 0
    dist_numeric <- 0
  }
  else {
    num <- length(idnum)
    msd <- mean(apply(data[, idnum, drop = FALSE], 2, sd))
    x <- as.matrix(data[, idnum, drop = FALSE])
    x <- scale(x, center = TRUE, scale = TRUE)
    dist_numeric <- switch(method,
                           gower = distNumeric(x, x, method = "mrw"),
                           wishart = distNumeric(x, x,  method = "sev"),
                           podani = distNumeric(x, x, method = "ser.2"),
                           huang = distNumeric(x, x, method = "se"),
                           harikumar = as.matrix(dist(x, method = "manhattan")),
                           ahmad = distNumeric(x, x, method = "se"))
  }
  if (is.null(idbin)) {
    bin <- 0
    dist_binary <- 0
  }
  else {
    bin <- length(idbin)
    dist_matchbin <- matching(data[, idbin, drop = FALSE],
                              data[, idbin, drop = FALSE])
    if (method == "ahmad") {
      dist_binary <- cooccur(data[, c(idbin, idcat), drop = FALSE])
    }
    else {
      if (method == "huang" | method == "harikumar") {
        dist_binary <- dist_matchbin * bin
      }
      else {
        dist_binary <- dist_matchbin
      }
    }
  }
  if (is.null(idcat)) {
    cat <- 0
    dist_cat <- 0
  }
  else {
    cat <- length(idcat)
    dist_matchcat <- matching(data[, idcat, drop = FALSE],
                              data[, idcat, drop = FALSE])
    if (method == "harikumar") {
      dist_cat <- cooccur(data[, idcat, drop = FALSE])
    }
    else {
      if (method == "huang") {
        dist_cat <- dist_matchcat * cat
      }
      else {
        if (method == "ahmad") {
          dist_cat <- dist_binary
        }
        else {
          dist_cat <- dist_matchcat
        }
      }
    }
  }
  nvar <- num + bin + cat
  dist_mix <- switch(method,
                     gower = dist_numeric * 1/nvar + dist_binary * bin/nvar + dist_cat * cat/nvar,
                     wishart = (dist_numeric * 1/nvar + dist_binary * bin/nvar + dist_cat * cat/nvar)^0.5,
                     podani = (dist_numeric + dist_binary * bin + dist_cat *
                                 cat)^0.5,
                     huang = dist_numeric + dist_binary * msd +
                       dist_cat * msd,
                     harikumar = dist_numeric + dist_binary +
                       dist_cat,
                     ahmad = dist_numeric + (dist_binary)^2)
  return(dist_mix)
}



#' Phyloregions for functional traits
#'
#' Generates a sparse community matrix as input for clustering regions
#' based on the similairity of functional traits across species.
#' @param x A community data in long format with one column representing
#' sites labeled \dQuote{grids} and another column representing species
#' labeled \dQuote{species}.
#' @param trait A data frame or matrix object with the first column
#' labeled \dQuote{species} containing the texonomic groups to be
#' evaluated whereas the remaining columns have the various functional
#' traits. The variables must be a mix of numeric and categorical values.
#' @param num A vector of column index of the numerical trait variable.
#' @param bin A vector of column index of the binary trait variable.
#' @param cat A vector of column index of the categorical trait variable.
#' @param cut The slice time.
#' @param phy is a dated phylogenetic tree with branch lengths stored
#' as a phylo object (as in the ape package).
#' @param k The desired number of clusters for single numeric trait
#' variables.
#' @param method Whether to compute phylogenetic
#' (method = \dQuote{phylo}), functional (method = \dQuote{trait}),
#' single functional categorical variables (method = \dQuote{single_cat}),
#' or for single numeric variables (method = \dQuote{single_num}).
#' @rdname trait_matrix
#' @keywords bioregion
#' @importFrom dbscan dbscan
#' @importFrom stats as.dist var sd dist
#' @importFrom utils stack
#' @importFrom Matrix Matrix
#' @importFrom ape keep.tip
#' @return Function returns a community data frame that captures the
#' count of each species based on its cluster membership.
#' @export
trait_matrix <- function (x, trait, num = NULL, bin = NULL,
                          cat = NULL, cut = NULL, phy = NULL,
                          method="trait", k = NULL)
{

  x$species <- gsub(" ", "_", x$species)

  if (!is.null(phy)) {
    subphy <- keep.tip(phy, intersect(phy$tip.label, x$species))
    submat <- subset(x, x$species %in% intersect(phy$tip.label, x$species))
    tx <- get_clades(subphy, cut = cut)
    z <- length(tx)
    memb <- rep(seq_len(z), lengths(tx))
    sp <- sparseMatrix(seq_along(memb), j = memb, dims=c(length(memb),z),
                       dimnames=list(unlist(tx), as.character(seq_len(z))))
  } else if(method=="trait") {
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    index <- intersect(x$species, trait$species)
    trait <- subset(trait, trait$species %in% index)

    row.names(trait) <- trait$species
    g <- gower_distance(trait, method = "gower",
                        idnum=num, idbin=bin, idcat=cat)
    gc()
    g1 <- as.dist(g)
    g2 <- dbscan(g1, minPts = 5, eps = 0.05)
    gc()

    memb <- g2$cluster
#    names(memb) <- labels(g1)
#    memb <- memb[!(memb==0)]
    z <- length(unique(memb))

    sp <- sparseMatrix(seq_along(memb), j = memb, dims=c(length(memb),z),
                       dimnames=list(labels(g1), as.character(seq_len(z))))

    submat <- subset(x, x$species %in% index)
  } else if (method=="single_cat"){
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    index <- intersect(x$species, trait$species)
    trait <- subset(trait, trait$species %in% index)

    row.names(trait) <- trait$species
    zz <- trait[, cat, drop=FALSE]
    names(zz) <- "foo"

    zz$foo[zz$foo==""] <- "noise"
    memb <- as.integer(factor(zz[,1]))
    z <- length(unique(memb))
    sp <- sparseMatrix(seq_along(memb), j = memb, dims=c(length(memb),z),
                       dimnames=list(row.names(zz), as.character(seq_len(z))))

    submat <- subset(x, x$species %in% index)
  } else if (method=="single_num"){
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    index <- intersect(x$species, trait$species)
    trait <- subset(trait, trait$species %in% index)

    row.names(trait) <- trait$species
    zz <- trait[, num, drop=FALSE]
    names(zz) <- "foo"
    zz$cluster <- findInterval(zz$foo,
                               unique(quantile(zz$foo, probs = seq(0, 1, 1/k))),
                               rightmost.closed = TRUE)
    memb <- zz$cluster
    z <- length(unique(memb))
    sp <- sparseMatrix(seq_along(memb), j = memb, dims=c(length(memb),z),
                       dimnames=list(row.names(zz), as.character(seq_len(z))))

    submat <- subset(x, x$species %in% index)
  }

  M <- long2sparse(submat)

  mx <- M %*% sp
  raw_mtx <- mx/rowSums(mx)
  res <- sparse2long(mx)
  return(list(comm_dat = res, k = z, raw_mtx = raw_mtx))
}
