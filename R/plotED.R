legend <- function(y, cols, vals, lab="ED", leg=5, lwd=15, pos="bottomright") {
  if(pos=="bottomright"){
    a=y@bbox[3]+0.5
    b=y@bbox[2]
  }
  if(pos=="topleft"){
    a=y@bbox[1]-0.5
    b=y@bbox[4]-leg
  }
  if(pos=="bottomleft"){
    a=y@bbox[1]-0.5
    b=y@bbox[2]
  }
  if(pos=="topright"){
    a=y@bbox[3]+0.5
    b=y@bbox[4]-leg
  }
  X <- cbind(rep(a, length(cols)), rep(a, length(cols)))
  Y <- b + cbind(0:(length(cols) - 1)/length(cols),
                 1:length(cols)/length(cols)) * (leg)
  for (i in 1:length(cols)) lines(X[i, ], Y[i, ], col = cols[i],
                                  lwd = lwd, lend = 2)
  text(x = a, y = b, round(min(vals), 3), pos = 4,
       cex = 0.7) # lim texts
  text(x = a, y = b+(leg/2), round(median(vals), 3), pos = 4,
       cex = 0.7)
  text(x = a, y = b+leg, round(max(vals), 3),
       pos = 4, cex = 0.7)
  text(x = a, y = b+leg, lab,
       pos = 3, cex = 1)
}

#' Map evolutionary distinctiveness of phyloregions in space
#'
#' This function maps evolutionary distinctiveness of bioregions based on the mean
#' value of phylogenetic beta diversity between a focal bioregion and all
#' other bioregions in the study area.
#' @param dat a pairwise matrix of beta diversity (phylo or non-phylo).
#' @param shp a polygon shapefile of grid cells.
#' @param method the agglomeration method to be used. This can be one
#' of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
#' "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param cex character (or symbol) expansion: a numerical vector.
#' @param k the optimal number of clusters derived from the elbow method
#' (as in the \code{GMD} package).
#' @param border color of the border.
#' @param swatch name of the palette to generate colors from. The name is matched
#' to the list of available color palettes from the code{hcl.colors} function in
#' the code{grDevices} package.
#' @param legend logical indicating whether to add a legend to the map.
#' @param pos location to position the legend such as \dQuote{bottomright},
#' \dQuote{bottomleft}, \dQuote{topleft}, and \dQuote{topright}.
#' @param title a character string indicating the caption to be placed on top
#' of the legend.
#' @param \dots arguments passed among methods.
#' @rdname plot_evoldistinct
#' @return Returns no value, just plot the evolutionary distinctiveness of bioregions
#' @importFrom stats hclust cutree
#' @importFrom rgeos gUnaryUnion
#' @importFrom sp SpatialPolygonsDataFrame merge
#' @importFrom sp CRS proj4string
#' @importFrom raster text
#' @importFrom graphics legend par points rect segments strheight strwidth text xinch yinch plot lines
#' @importFrom grDevices rgb hcl.colors as.graphicsAnnot xy.coords
## @method plot default
#'
#' @inheritParams evoldistinct
#' @inheritParams choropleth

#' @author Barnabas H. Daru

#' @references
#' Daru, B.H., Van der Bank, M., Maurin, O., Yessoufou, K., Schaefer, H., Slingsby, J.A. & Davies, T.J. (2016) A novel phylogenetic
#' regionalization of the phytogeographic zones of southern Africa reveals their hidden evolutionary affinities.
#' \emph{Journal of Biogeography} \strong{43}: 155-166.
#'
#' Daru, B.H., Elliott, T.L., Park, D.S. & Davies, T.J. (2017) Understanding the processes underpinning patterns of phylogenetic
#' regionalization. \emph{Trends in Ecology and Evolution} \strong{32}: 845-860.
#' @examples
#' x <- c("phyloregion", "raster", "Matrix", "ape", "betapart",
#' "rgeos", "vegan", "colorspace")
#' lapply(x, require, character.only = TRUE) # load the required packages
#'
#' data(africa)
#' tree <- africa$phylo
#' x <- sampl2sparse(africa$comm)
#'
#' subphy <- match_phylo_comm(tree, x)$phy
#' submat <- match_phylo_comm(tree, x)$com
#'
#' pbc <- phylobeta(submat, subphy)
#' plot_evoldistinct(pbc, shp=africa$polys)
#' @export plot_evoldistinct
plot_evoldistinct <- function (dat,
                               method = "average",
                               shp = shp,
                               k = 10,
                               cex = 1,
                               swatch = "YlOrBr",
                               pos = "bottomright",
                               lab = "",
                               legend = TRUE,
                               border = NA,
                               leg =5,
                               lwd = 15, ...) {
  if (!inherits(dat, "matrix")) {
    stop("dat must be of class 'matrix'")
  }
  P <- dat
  Q = as.dist(P)
  P1 <- hclust(Q, method = method)
  g <- cutree(P1, k)
  dx <- data.frame(cluster = g)
  dx <- data.frame(grids = row.names(dx), dx)
  ED <- evoldistinct(P, k, method = method)
  m <- sp::merge(shp, dx, by = "grids")
  if (!inherits(m, "SpatialPolygons")) {
    stop("Invalid geometry, may only be applied to polygons")
  }
  m <- m[!is.na(m@data$cluster), ]
  region <- rgeos::gUnaryUnion(m, id = m@data$cluster)
  row.names(region) <- as.character(1:length(region))
  fx <- unique(m$cluster)
  fx <- as.data.frame(fx)
  colnames(fx) <- "cluster"
  region <- SpatialPolygonsDataFrame(region, fx)
  names(ED[[1]])[2] <- "cluster"
  m1 <- sp::merge(region, ED[[1]], by = "cluster")
  proj4string(m1) = proj4string(shp)
  COLRS <- rev(hcl.colors(k, swatch))
  y <- choropleth(m1, m1$ED, k) #, style = style
  plot(y, col = COLRS[y$values], border = border, ...)
  text(y, labels = as.character(y$cluster), cex = cex)
  if(legend){
    legend(y,COLRS,vals = m1$ED,leg = leg,lwd = lwd,pos = pos,lab = lab)
  }
}
