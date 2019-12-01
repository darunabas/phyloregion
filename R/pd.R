#' Phylogenetic diversity
#'
#' \code{PD.sparse} calculates Faith's phylogenetic diversity and much faster
#' using sparse matrix.
#'
#' @param x a community matrix, i.e. an object of class matrix or Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @return a vector with the PD for all samples.
#' @keywords cluster
#' @seealso read.community read.tree phylobeta_core
#' @examples
#' library(ape)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), tree$tip.label))
#' PD.sparse(com, tree)
#' @rdname PD.sparse
#' @export
PD.sparse <- function(x, phy){
  y <- phylo_community(x, phy)
  el <- attr(y, "edge.length")
  res <- vapply(y, function(x, el)sum(el[x]), 0, el)
  res
}

