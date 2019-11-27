#' Plot phyloregions in NMDS ordinatoion space
#' 
#' This function plots clusters of cells in NMDS space.
#' @param x a pairwise matrix of beta diversity (phylo or non-phylo).

#' @param method the agglomeration method to be used. This should 
#' be (an unambiguous abbreviation of) one of \dQuote{ward.D}, \dQuote{ward.D2}, 
#' \dQuote{single}, \dQuote{complete}, \dQuote{average} (= UPGMA), \dQuote{mcquitty} (= WPGMA), 
#' \dQuote{median} (= WPGMC) or \dQuote{centroid} (= UPGMC).
#' @param k the optimal number of clusters derived from the elbow method 
#' (as in the \code{GMD} package).
#' @param ... graphical parameters to plot
#' @rdname plot.NMDS
#' @return Returns no value, just plot the bioregions
#' @importFrom colorspace hex2RGB
#' @importFrom vegan metaMDS
#' @importFrom graphics plot Axis
#' @method plot default
#' @keywords visualization and mapping
#' 
#' @inheritParams evoldistinct 

#' @author Barnabas H. Daru

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
#' plot.NMDS(pbc)
#' @export plot.NMDS

plot.NMDS <- function(x, 
                      k = 10, 
                      method = "average", ...){
  
  ED <- evoldistinct(x, 
                     k=k, 
                     method = method)
  c1 <- vegan::metaMDS(ED[[2]], trace = 0)
  

  plot(c1$points, 
       pch = 21, 
       cex = 8, 
       bg = hexcols(c1), 
       lwd=2, 
       xlab="NMDS axis 1", 
       ylab="", 
       axes=FALSE, 
       frame.plot=FALSE, 
       las=1, 
       xlim=c(min(c1$points-0.05), max(c1$points+0.05)), 
       ylim=c(min(c1$points+0.05), max(c1$points+0.05)), 
       yaxs="r", 
       xaxs="r", 
       cex.lab=1.5, ...)
  
  box(which = "plot", bty = "l", lwd=2)
  
  Axis(side=1, labels=TRUE, lwd=2, las=1, cex.axis=1.2)
  Axis(side=2, labels=TRUE, lwd=2, las=1, cex.axis=1.2)
  
  text(c1$points, as.character(1:c1$nobj), cex =2)
  
}




