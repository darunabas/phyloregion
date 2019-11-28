#' Slice phylogenetic tree at various time depths
#'
#' This function slices a dated phylogenetic tree at successive time depths back
#' in time by collapsing younger phylogenetic branches into older ones to infer
#' the origins of species assemblages.
#'
#' @param phy A dated phylogenetic tree as an object of class \dQuote{phylo}.
#' @param n Time depth to slice the phylogenetic tree (often in millions of years
#' for dated trees).
#' @param collapse Logical, collapse internal edges with zero edge length.
#' @param \dots arguments passed among methods.
#' @rdname timeslice
#' @keywords phyloregion
#' @importFrom ape node.depth.edgelength di2multi
#'
#' @return
#' {A tree with the phylogenetic structure removed at the specified time depth}
#'
#' @references
#' Daru, B.H., van der Bank, M. & Davies, T.J. (2018) Unravelling the evolutionary origins of biogeographic
#' assemblages. \emph{Diversity and Distributions} \strong{24}: 313–324.
#'
#'
#' @author Barnabas H. Daru \email{darunabas@gmail.com}
#'
#' @examples
#' require(ape)
#'
#' set.seed(1)
#' tree <- rcoal(50)
#' x <- timeslice(tree,.5)
#'
#' par(mfrow=c(1,2))
#' plot(tree); axisPhylo()
#' plot(x); axisPhylo()
#' @export
timeslice <- function(phy, n=0.2, collapse=FALSE, ...){
    node.age <- node.depth.edgelength(phy)
    max.age <- max(node.age)
    if(n > max.age) stop("n is too large!")
    goal.length <- max.age - n
    phy$edge.length[node.age[phy$edge[,1]] > goal.length] <- 0
    ind <- which( (node.age[phy$edge[,1]] < goal.length) &
                  (node.age[phy$edge[,2]] > goal.length) )
    phy$edge.length[ind] <- goal.length - node.age[phy$edge[ind,1]]
    if(collapse) phy <- di2multi(phy)
    phy
}

