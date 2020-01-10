#' Get descendant nodes of phylogeny at a given time depth
#'
#' \code{get_clades} returns the tips that descend from a given node or time depth on a dated phylogenetic tree.
#'
#' @param tree is a dated phylogenetic tree with branch lengths stored
#' as a phylo object (as in the \code{ape} package).
#' @param ... Further arguments passed to or from other methods.
#' @param cut the slice time
#' @rdname get_clades
#' @keywords bioregion
#' @importFrom phangorn Descendants
#' @importFrom ape node.depth.edgelength
#'
#' @return A list of descendants
#'
#' @references
#' Schliep, K.P. (2010) phangorn: phylogenetic analysis in
#' R. \emph{Bioinformatics} \strong{27}: 592–593.
#'
#' @examples
#' require(ape)
#' data(bird.orders)
#' plot(bird.orders)
#' axisPhylo(side = 1)
#' abline(v=28-23) # the root is here at 28
#' get_clades(bird.orders, 23)
#' @export
get_clades <- function(tree, cut=2, ...){
  nh <- node.depth.edgelength(tree)
  nh <- max(nh) - nh
  ind <- which( (nh[tree$edge[,1]] > cut) & (nh[tree$edge[,2]] <= cut) )
  desc <- Descendants(tree)
  res <- desc[tree$edge[ind,2]]
  lapply(res, function(res, tips)tips[res], tree$tip.label)
}

