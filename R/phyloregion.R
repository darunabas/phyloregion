dissolve_poly <- function(x){
  # Now the dissolve
  x <- x[!is.na(x@data$cluster),]
  region <- rgeos::gUnaryUnion(x, id = x@data$cluster)
  # make sure row names match
  row.names(region) <- as.character(seq_along(region))
  # Extract the data you want (the larger geography)
  fx <- unique(x$cluster)
  fx <- as.data.frame(fx)
  colnames(fx) <- "cluster"
  # And add the data back in
  SpatialPolygonsDataFrame(region, fx)
}
#' Calculate evolutionary distinctiveness of phyloregions
#'
#' This function estimates evolutionary distinctiveness of each phyloregion by
#' computing the mean value of phylogenetic beta diversity between a focal
#' phyloregion and all other phyloregions in the study area.
#'
#' @param x A distance matrix
#' @param k The desired number of phyloregions, often as determined by
#' \code{optimal_phyloregion}.
#' @param method the agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of \dQuote{ward.D}, \dQuote{ward.D2}, \dQuote{single},
#' \dQuote{complete}, \dQuote{average} (= UPGMA), \dQuote{mcquitty} (= WPGMA), \dQuote{median}
#' (= WPGMC) or \dQuote{centroid} (= UPGMC).
#' @param shp a polygon shapefile of grid cells.
#' @param ... Further arguments passed to or from other methods.
#' @rdname phyloregion
#' @keywords phyloregion
#' @importFrom stats as.dist hclust cutree
#' @importFrom rgeos gUnaryUnion
#' @importFrom sp SpatialPolygonsDataFrame merge
#' @importFrom sp CRS proj4string
#' @importFrom raster text
#' @importFrom graphics legend par points rect segments strheight strwidth text
#' xinch yinch lines
#' @importFrom grDevices rgb hcl.colors as.graphicsAnnot xy.coords
#'
#' @return
#' An object of class \code{phyloregion} containing
#' \itemize{
#'   \item a data frame membership with columns grids and cluster
#'   \item k the number of clusters
#'   and additionally there can be an shape file and other bjects.
#'   This representation may still change.
#' }
#' @references
#'
#' Daru, B.H., Van der Bank, M., Maurin, O., Yessoufou, K., Schaefer, H.,
#' Slingsby, J.A. & Davies, T.J. (2016) A novel phylogenetic regionalization of
#' the phytogeographic zones of southern Africa reveals their hidden
#' evolutionary affinities. \emph{Journal of Biogeography} \strong{43}: 155-166.
#'
#' Daru, B.H., Elliott, T.L., Park, D.S. & Davies, T.J. (2017) Understanding the
#' processes underpinning patterns of phylogenetic regionalization.
#' \emph{Trends in Ecology and Evolution} \strong{32}: 845-860.
#'
#' Daru, B.H., Holt, B.G., Lessard, J.P., Yessoufou, K. & Davies, T.J. (2017)
#' Phylogenetic regionalization of marine plants reveals close evolutionary
#' affinities among disjunct temperate assemblages.
#' \emph{Biological Conservation} \strong{213}: 351-356.

#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#' @seealso \code{\link{evol_distinct}}, \code{\link{optimal_phyloregion}},
#' \code{\link[picante]{evol.distinct}} for a different approach.
#'
#' @examples
#' library(ape)
#' tree <- read.tree(text = "((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
#'   c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
#'   dimnames = list(paste0("g", 1:6), tree$tip.label))
#' pbc <- phylobeta(com, tree)
#' phyloregion(pbc[[1]], k = 3)
#' @export
phyloregion <- function(x, k = 10, method = "average", shp = NULL, ...) {

  Q <- as.dist(x)
  P1 <- hclust(Q, method = method)
  g <- cutree(P1, k)
  dx <- data.frame(grids=names(g), cluster = unname(g))

  x <- as.matrix(x)
  colnames(x) <- rownames(x)

  region.mat <- matrix(NA, k, k, dimnames = list(1:k, 1:k))

  for (i in 1:k) {
    for (j in 1:k) {
      region.mat[i, j] <- mean(x[names(g)[g == i], names(g)[g == j]])
    }
  }
  region.dist <- as.dist(region.mat)
  region.mat <- as.matrix(region.dist)

  evol_distinct <- colSums(region.mat) / (nrow(region.mat) - 1)

  evol_distinct <- data.frame(ED = evol_distinct)
  evol_distinct <- cbind(cluster = rownames(evol_distinct),
    data.frame(evol_distinct, row.names = NULL))

  if (length(shp) == 0) {
    r <- list(membership=dx, k=k,
      evol_distinct = evol_distinct, region.dist = region.dist,
      region.df = dx)
    class(r) <- c("phyloregion")
    r
  } else {

    m <- sp::merge(shp, dx, by = "grids")
    if (!inherits(m, "SpatialPolygons")) {
      stop("Invalid geometry, may only be applied to polygons")
    }
    m <- m[!is.na(m@data$cluster), ]

    region <- dissolve_poly(m)

    m1 <- sp::merge(region, evol_distinct, by = "cluster")
    proj4string(m1) <- proj4string(shp)

    c1 <- vegan::metaMDS(region.dist, trace = 0)

    v <- data.frame(hex2RGB(hexcols(c1))@coords)
    v$r <- v$R * 255
    v$g <- v$G * 255
    v$b <- v$B * 255

    v$COLOURS <- rgb(v$r, v$g, v$b, maxColorValue = 255)
    v$cluster <- rownames(v)

    y <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE),
      list(region, v, m1))

    index <- match(dx$cluster, y$cluster)
    z <- cbind(dx, ED = y$ED[index], COLOURS = y$COLOURS[index])

    # membership

    r <- list(membership=dx, k=k, shp = y,
              region.dist = region.dist, region.df = z, NMDS = c1)
    class(r) <- "phyloregion"
    r
  }
}


#' @rdname phyloregion
#' @importFrom igraph graph_from_incidence_matrix cluster_infomap communities
#' @importFrom igraph membership
#' @export
infomap <- function(x, shp = NULL, ...){
  x@x[x@x >1e-8] <- 1
  g <- graph_from_incidence_matrix(x)
  imc <- cluster_infomap(g, ...)
  ms <- membership(imc)
  k <- max(ms)
  ind <- names(ms) %in% rownames(x)
  dx <- data.frame(grids=names(ms)[ind], cluster=unname(ms)[ind])
  if(!is.null(shp)){
     shp <- sp::merge(shp, dx, by = "grids")
     shp <- dissolve_poly(shp)
  }
  result <- list(membership=dx, k=k, shp=shp)
  class(result) <- "phyloregion"
  result
}
