#' A function to calculate evolutionary distinctiveness of bioregions.
#'
#' This function estimates evolutionary distinctiveness of each bioregion by
#' computing the mean value of phylogenetic beta diversity between a focal
#' bioregion and all other bioregions in the study area.
#'
#' @param x A distance matrix
#' @param k The desired number of bioregions
#' @param method the agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#' "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC).
#' @param ... Further arguments passed to or from other methods.
#' @rdname evoldistinct
#' @keywords phyloregion
#' @importFrom stats as.dist hclust cutree
#'
#' @return
#' \item{psim}{A bioregion × bioregion phylogenetic beta diversity distance matrix}
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
#' @seealso \code{\link[picante]{evol.distinct}} for a different approach.
#'
#' @examples
#' library(ape)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), tree$tip.label))
#' pbc <- phylobeta(com, tree)
#' evoldistinct(pbc, k=3)
#' @export
evoldistinct <- function(x, k=10, method="average", ...){
  P1 <- as.dist(x)

  P2 <- hclust(P1, method="average")
  g <- cutree(P2, k)

  d <- data.frame(bioreg=g)
  d$grids <- rownames(d)

  x <- as.matrix(x)
  colnames(x) <- rownames(x)

  region.mat <- matrix(NA, k, k, dimnames = list(1:k, 1:k))

  for(i in 1:k){

    for(j in 1:k){
      region.mat[i, j] <- mean(x[names(g)[g == i], names(g)[g == j]])
    }

  }
  region.dist <- as.dist(region.mat)
  region.mat <- as.matrix(region.dist)

  evol_distinct <- colSums(region.mat)/(nrow(region.mat)-1)

  evol_distinct <- data.frame(ED=evol_distinct)
  evol_distinct$cluster <- rownames(evol_distinct)
  return(list(evol_distinct=evol_distinct, region.dist=region.dist))
}
