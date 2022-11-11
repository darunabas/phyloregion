#' Get descendant nodes of phylogeny at a given time depth
#'
#' \code{get_clades} returns the tips that descend from a given node or time
#' depth on a dated phylogenetic tree.
#'
#' @param tree is a dated phylogenetic tree with branch lengths stored
#' as a phylo object (as in the \code{ape} package).
#' @param cut the slice time
#' @param k number of slices
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
get_clades <- function(tree, cut=NULL, k=NULL){
  nh <- node.depth.edgelength(tree)
  nh <- max(nh) - nh
  if(!is.null(k)){
    if(k >= Ntip(tree)) return(as.list(tree$tip.label))
    if(k == 1) return(list(tree$tip.label))
    kids <- lengths(Descendants(tree, type = "children"))
    kids[kids>0] <- kids[kids>0]- 1L
    tmp <- 1
    eps <- 1e-8
    ordered_nh <- order(nh, decreasing = TRUE)
    i <- 1
    while(tmp < k){
      j <- ordered_nh[i]
      cut <- nh[j] - eps
      tmp <- tmp + kids[j]
      i <- i+1
    }
  }
  ind <- which( (nh[tree$edge[,1]] > cut) & (nh[tree$edge[,2]] <= cut) )
  desc <- Descendants(tree)
  res <- desc[tree$edge[ind,2]]
  lapply(res, function(res, tips)tips[res], tree$tip.label)
}

