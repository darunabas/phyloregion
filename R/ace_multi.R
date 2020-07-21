#' Ancestral Character Estimation for Multiple Traits
#'
#' This function reconstructs ancestral states for multiple
#' continuous (numeric) trait using phylogenetic independent
#' contrasts (PIC; Felsenstein, 1985).
#'
#' @param phy A dated and rooted phylogenetic tree of the class "phylo".
#' The root is assumed to be the unique node with no incoming edge.
#' @param x A species \code{×} k matrix, where the column sums equal 1.
#' These values represent the known states or contributions of each
#' species to each cluster/biome to be represented in the tree.
#' @importFrom castor asr_independent_contrasts
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{ancestral_states} The estimates of the ancestral
#'   character values.
#'   \item \code{CI95} The estimated 95 percent confidence intervals.
#'   \item \code{se} The standard-errors of estimated rates.
#'   \item \code{tree} The pruned tree with tips matching the input
#'   matrix \code{x} above.
#'   \item \code{piecolors} List of HUE colors for plotting nodepies in
#'   the tree.
#' }
#' @references Felsenstein, J. (1985) Phylogenies and the comparative
#' method. \emph{American Naturalist} \strong{125}: 1-15.
#' @examples
#' \donttest{
#' old.par <- par(no.readonly = TRUE)
#' require(ape)
#' data(africa)
#' m <- ace_multi(africa$theta, africa$phylo)
#'
#' par(mai=rep(0,4))
#' par(mfrow=c(2,1))
#' plot(m$tree, type = "fan", show.tip.label=FALSE,
#'      open.angle = 180, edge.width=0.5)
#'
#' K <- ncol(africa$omega)
#' COLRS <- phyloregion:::hue(K)
#'
#' nodepie(m$ancestral_states, radius = 2, col=COLRS)
#'
#'
#' plot_pie(africa$omega, shp = africa$polys, col=COLRS)
#' }
#'
#' @export
ace_multi <- function(x, phy) {

  row.names(x) <- gsub(" ", "_", row.names(x))
  index <- intersect(row.names(x), phy$tip.label)
  x <- x[index,]

  tree <- keep.tip(phy, index)

  M1 <- apply(x, 2, function(i) {
    i <- as.numeric(i)
    y <- asr_independent_contrasts(tree, i,
                                   weighted=TRUE,
                                   include_CI=TRUE)$ancestral_states
    y/sum(y)
  })
  M2 <- apply(x, 2, function(i) {
    i <- as.numeric(i)
    y <- asr_independent_contrasts(tree, i,
                                   weighted=TRUE,
                                   include_CI=TRUE)$CI95
    y/sum(y)
  })
  M3 <- apply(x, 2, function(i) {
    i <- as.numeric(i)
    y <- asr_independent_contrasts(tree, i,
                                   weighted=TRUE,
                                   include_CI=TRUE)$standard_errors
    y/sum(y)
  })
  COLRS <- hue(ncol(x))
  return(list(ancestral_states = M1, CI95 = M2, se = M3,
              tree = tree, piecolors = COLRS))
}

