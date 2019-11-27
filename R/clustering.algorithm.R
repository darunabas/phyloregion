#' Cluster algorithm selection and validation
#'
#' This function contrasts different hierarchical clustering algorithms 
#' on the phylogenetic beta diversity matrix for degree of 
#' data distortion using Sokal & Rohlf’s (1962) cophenetic 
#' correlation coefficient.
#'
#' @param x a numeric matrix, data frame or "dist" object.
#' @param rot_angle rotation angle of the barplot x-axis labels.
#' @param cex character (or symbol) expansion: a numerical vector for 
#' the size of the x-axis labels
#' @param \dots arguments passed among methods.
#' @rdname clustering.algorithm
#' @keywords bioregion
#' @importFrom stats hclust as.dist cor cophenetic
#' @importFrom cluster diana
#' @importFrom graphics barplot
#' 
#' @references 
#' Sokal, R.R. & Rohlf, F.J. (1962) The comparison of dendrograms by 
#' objective methods. \emph{Taxon} \strong{11}: 33–40.
#'
#' @return
#' \itemize{
#'   \item A numeric value corresponding to the good clustering algorithm 
#'   for the distance matrix
#'   \item If plot = TRUE, a barplot of cophenetic correlation for all
#'   the clustering algorithms is drawn.
#' }
#'
#' @examples
#' require(betapart)
#' data(africa)
#' tree <- africa$phylo
#' M <- sampl2sparse(africa$comm, method = "nonphylo")
#' 
#' bc <- beta.core(M)
#' clustering.algorithm(bc, plot = TRUE)
#' @export clustering.algorithm
clustering.algorithm <- function(x, plot=FALSE, 
                                 rot_angle=45, cex=0.7, col="grey", ...){
  hc1 <- hclust(as.dist(x), method="average")
  hc2 <- hclust(as.dist(x), method="single")
  hc3 <- hclust(as.dist(x), method="complete")
  hc4 <- hclust(as.dist(x), method="ward.D")
  hc5 <- hclust(as.dist(x), method="ward.D2")
  hc6 <- hclust(as.dist(x), method="mcquitty")
  hc7 <- hclust(as.dist(x), method="median")
  hc8 <- hclust(as.dist(x), method="centroid")
  hc9 <- diana(as.dist(x))
  

  z <- c(UPGMA=cor(as.dist(x), cophenetic(hc1)), 
    Single=cor(as.dist(x), cophenetic(hc2)),
    Complete=cor(as.dist(x), cophenetic(hc3)),
    ward.D=cor(as.dist(x), cophenetic(hc4)),
    ward.D2=cor(as.dist(x), cophenetic(hc5)),
    WPGMA=cor(as.dist(x), cophenetic(hc6)),
    WPGMC=cor(as.dist(x), cophenetic(hc7)),
    UPGMC=cor(as.dist(x), cophenetic(hc8)),
    DIANA=cor(as.dist(x), cophenetic(hc9)))
  plt <- barplot(z, las=1, 
                 ylab = "Cophenetic correlation", 
                 plot = plot, xaxt="n", col = col,
                 ylim=c(0,0.1+max(z)), ...)
  text(plt, par("usr")[3], 
       labels = names(z), 
       srt = rot_angle, 
       adj = c(1.1,1.1), 
       xpd = TRUE, cex=cex) 
  
  cat("\nA good clustering algorithm for the distance matrix is:\n", 
      names(z[which.max(z)]), 
      "with cophenetic correlation =", 
      z[which.max(z)], 
      "\n\n")
  
  return(z)

}

