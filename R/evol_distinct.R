out_degree <- function(x){
  res <- numeric(max(x$edge))
  res[seq_len(Ntip(x))] <- 1L
  for(i in seq_len(nrow(x$edge))) res[x$edge[i,1]] <- res[x$edge[i,1]] + 1L
  res
}

sons <- function(x){
  x <- reorder(x, 'postorder')
  res <- numeric(max(x$edge))
  res[seq_len(Ntip(x))] <- 1L
  for(i in seq_len(nrow(x$edge))) {
    res[x$edge[i,1]] <- res[x$edge[i,1]] + res[x$edge[i,2]]
  }
  res
}


#' Species' evolutionary distinctiveness
#'
#' Calculates evolutionary distinctiveness measures for a suite of species by:
#' a) equal splits (Redding and Mooers 2006)
#' b) fair proportions (Isaac et al., 2007).
#' This a new implementation of the picante function \code{evol.distinct}
#' however allowing multifurcations and can be orders of magnitude faster.
#'
#' @param tree an object of class \code{phylo}.
#' @param type a) equal splits (Redding and Mooers 2006) or b) fair proportions
#' (Isaac et al., 2007)
#' @param scale The scale option refers to whether or not the phylogeny should
#' be scaled to a depth of 1 or, in the case of an ultrametric tree, scaled such
#' that branch lengths are relative.
#' @param use.branch.lengths If use.branch.lengths=FALSE, then all branch
#' lengths are changed to 1.
#' @param \dots Further arguments passed to or from other methods.
#'
#' @seealso \code{\link[picante]{evol.distinct}}, \code{\link{phyloregion}}
#' @return a named vector with species scores.
#' @author Klaus Schliep
#' @references
#' Redding, D.W. and Mooers, A.O. (2006). Incorporating evolutionary measures
#' into conservation prioritisation. \emph{Conservation Biology}, \strong{20},
#' 1670--1678.
#'
#' Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M.
#' (2007). Mammals on the EDGE: conservation priorities based on threat and
#' phylogeny. \emph{PLoS ONE}, \strong{2}, e296.
#' @examples
#' tree <- ape::rcoal(10)
#' evol_distinct(tree)
#' evol_distinct(tree, type = "fair.proportion")
#' @importFrom ape is.rooted is.ultrametric reorder.phylo branching.times
#' @importFrom stats reorder
#' @export
evol_distinct <- function (tree, type = c("equal.splits", "fair.proportion"),
    scale = FALSE, use.branch.lengths = TRUE, ...) {
    type <- match.arg(type)
    if (is.rooted(tree) == FALSE)
        warning("A rooted phylogeny is required for meaningful output of this function",
            call. = FALSE)
    if (scale == TRUE) {
        if (is.ultrametric(tree) == TRUE)
            tree$edge.length <- tree$edge.length /
              as.numeric(branching.times(tree)[1])
        else tree$edge.length <- tree$edge.length/sum(tree$edge.length)
    }
    if (use.branch.lengths == FALSE)
        tree$edge.length <- rep(1, length(tree$edge.length))

    el <- res <- numeric(max(tree$edge))
    el[tree$edge[,2]] <- tree$edge.length

    tree <- reorder(tree)

    if(type == "equal.splits"){
      od <- out_degree(tree)
      for(i in seq_len(nrow(tree$edge))){
        res[tree$edge[i,2]] <- (res[tree$edge[i,1]] + el[tree$edge[i,2]]) /
            od[tree$edge[i,2]]
      }
    } else if (type == "fair.proportion") {
        so <- sons(tree)
        for(i in seq_len(nrow(tree$edge))){
            res[tree$edge[i,2]] <- res[tree$edge[i,1]] + el[tree$edge[i,2]] /
                so[tree$edge[i,2]]
        }
    } else return(NULL)
    res <- res[seq_len(Ntip(tree))]
    names(res) <- tree$tip.label
    return(res)
}
