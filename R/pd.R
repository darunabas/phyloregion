#' Phylogenetic diversity and faster
#'
#' \code{PD} calculates Faith's (1992) phylogenetic diversity and much faster
#' using sparse matrix.
#'
#' @param x a community matrix, i.e. an object of class matrix or Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @param method If \code{method =} \dQuote{comm} (the default),
#' the phylogenetic diversity of each community in the study area;
#' else \dQuote{total}, the total PD of the study area.
#' @return a vector with the PD for all samples.
#' @keywords cluster
#' @seealso read.community read.tree phylobeta_core
#' @references
#' Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity.
#' \emph{Biological Conservation} \strong{61}: 1–10.
#' @examples
#' library(ape)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), tree$tip.label))
#' PD(com, tree)
#' @rdname PD
#' @export
PD <- function(x, phy, method = "comm"){
    if (method == "comm") {
        y <- phylo_community(x, phy)
        el <- attr(y, "edge.length")
        res <- vapply(y, function(x, el)sum(el[x]), 0, el)
    }
    else if (method == "total") {
        res <- sum(phy$edge.length)
    }
  res
}

