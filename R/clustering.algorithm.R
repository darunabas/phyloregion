#' Cluster algorithm selection and validation
#'
#' This function contrasts different hierarchical clustering algorithms
#' on the phylogenetic beta diversity matrix for degree of
#' data distortion using Sokal & Rohlf’s (1962) cophenetic
#' correlation coefficient.
#'
#' @param x a numeric matrix, data frame or "dist" object.
#' @rdname select_linkage
#' @keywords bioregion
#' @importFrom stats hclust as.dist cor cophenetic
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
#' data(africa)
#' tree <- africa$phylo
#' bc <- beta_diss(africa$comm)
#' y <- select_linkage(bc[[1]])
#' barplot(y, horiz = TRUE, las = 1)
#'
#' @export select_linkage
select_linkage <- function(x) {

  dist_x <- as.dist(x)

  methods_hc <- c("average", "single", "complete",
                  "ward.D", "ward.D2", "mcquitty",
                  "median", "centroid")

  z <- sapply(methods_hc, function(y){

    clust_methods <- hclust(dist_x, method = y)

    cor(dist_x, cophenetic(clust_methods), use="complete.obs")

  } )

#  z <- setNames(z, c('UPGMA', 'Single', 'Complete',
#                     'ward.D', 'ward.D2', 'WPGMA',
#                     'WPGMC', 'UPGMC'))
  names(z) <- c('UPGMA', 'Single', 'Complete', 'ward.D', 'ward.D2', 'WPGMA',
                'WPGMC', 'UPGMC')

  cat("\nA good clustering algorithm for the distance matrix is:\n",
      names(z[which.max(z)]),
      "with cophenetic correlation =",
      z[which.max(z)],
      "\n\n")

  return(z)
}
