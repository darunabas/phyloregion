#' Map phyloregions and bioregions
#' 
#' This function plots clusters of cells i.e \sQuote{phyloregions} or \sQuote{bioregions} in 
#' geographic space based on turnover of species or phylogenetic lineages.
#' @param dat a pairwise matrix of beta diversity (phylo or non-phylo).

#' @param shp a polygon shapefile of grid cells.

#' @param method the agglomeration method to be used. This should 
#' be (an unambiguous abbreviation of) one of \dQuote{ward.D}, \dQuote{ward.D2}, 
#' \dQuote{single}, \dQuote{complete}, \dQuote{average} (= UPGMA), \dQuote{mcquitty} (= WPGMA), 
#' \dQuote{median} (= WPGMC) or \dQuote{centroid} (= UPGMC).
#' @param cex character (or symbol) expansion: a numerical vector.

#' @param k the optimal number of clusters derived from the elbow method 
#' (as in the \code{GMD} package).
#' @param border color of the border.
#' @param ... graphical parameters to plot
#' @rdname plot.phyloregion
#' @return Returns no value, just plot the bioregions
#' @importFrom stats hclust cutree
#' @importFrom rgeos gUnaryUnion
#' @importFrom sp SpatialPolygonsDataFrame merge
#' @importFrom vegan metaMDS
#' @importFrom colorspace hex2RGB
#' @importFrom grDevices rgb
#' @importFrom sp CRS proj4string plot
#' @importFrom raster text
#' @importFrom graphics plot
#' @method plot default
#' @keywords visualization and mapping
#' 
#' @inheritParams evoldistinct 

#' @author Barnabas H. Daru

#' @references 
#' Daru, B.H., Van der Bank, M., Maurin, O., Yessoufou, K., Schaefer, H., Slingsby, J.A. & Davies, T.J. (2016) A novel phylogenetic 
#' regionalization of the phytogeographic zones of southern Africa reveals their hidden evolutionary affinities. 
#' \emph{Journal of Biogeography} \strong{43}: 155-166.
#' 
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
#' pbc <- phylobeta.core(submat, subphy)
#' plot.phyloregion(pbc, shp=africa$polys)
#' @export plot.phyloregion


plot.phyloregion <- function(dat, 
                           method="average", 
                           shp=shp, 
                           k=15, 
                           border=NA, 
                           cex=1, ...){
  
  # checks
  if(!inherits(dat,"matrix")){
    stop("dat must be of class 'matrix'")
  }
  P <- dat
  Q=as.dist(P)
  
  P1 <- hclust(Q, method=method)
  
  g <- cutree(P1, k)
  
  dx <- data.frame(cluster=g)
  
  dx <- data.frame(grids = row.names(dx), dx)
  
  ED <- evoldistinct(P, k, method = method)
  #signature(shp = "SpatialPolygons")
  
  m <- sp::merge(shp, dx, by="grids")
  if (!inherits(m, "SpatialPolygons")){
    stop("Invalid geometry, may only be applied to polygons")
  }
  m <- m[!is.na(m@data$cluster),]
  # Now the dissolve
  region <- rgeos::gUnaryUnion(m, id = m@data$cluster)
  
  # make sure row names match
  row.names(region) <- as.character(1:length(region))
  
  # Extract the data you want (the larger geography)
  fx <- unique(m$cluster)
  fx <- as.data.frame(fx)
  colnames(fx) <- "cluster"
  
  # And add the data back in
  region <- SpatialPolygonsDataFrame(region, fx)
  #proj4string(region) <- proj4string(shp)
  names(ED[[1]])[2] <- "cluster"
  
  m1 <- sp::merge(region, ED[[1]], by="cluster")
  proj4string(m1) = proj4string(shp) 
  
  c1 <- vegan::metaMDS(ED[[2]], trace = 0)
  
  v <- data.frame(hex2RGB(hexcols(c1))@coords)
  v$r <- v$R*255
  v$g <- v$G*255
  v$b <- v$B*255
  
  v$COLOURS <- rgb(v$r, v$g, v$b, maxColorValue = 255)
  
  v$cluster <- rownames(v)
  y <- merge(region, v, by="cluster")
  y <- y[, c("cluster", "COLOURS")]
  plot(y, col=y$COLOURS, border = border, ...)
  raster::text(y, labels=as.character(y@data$cluster), cex = cex)
  
}



