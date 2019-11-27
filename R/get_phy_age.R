#' Get age of nodes on phylogenetic tree
#'
#' This function computes the crown ages of nodes on a dated phylogenetic tree.
#'
#' @param phy A phylogenetic tree of the class \dQuote{phylo}.
#' @param n An index or indices indicating the node positions.
#' @rdname getNage
#' @keywords bioregion
#' @importFrom ape node.depth.edgelength
#' @return A numeric vector indicating the ages of the nodes.
#' @examples
#' require(ape)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' plot(tree)
#' nodelabels()
#' axisPhylo()
#' getNage(tree, 5)
#' @export
getNage <- function(phy, n){
  nd <-node.depth.edgelength(phy)
  nd <- max(nd) - nd
  nd[n]
}

