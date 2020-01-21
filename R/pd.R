#' Phylogenetic diversity
#'
#' \code{PD} calculates Faith's (1992) phylogenetic diversity.
#'
#' @param x a community matrix, i.e. an object of class matrix or Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @return a vector with the PD for all samples.
#' @keywords cluster
#' @seealso read.community read.tree phylobeta_core
#' @references
#' Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity.
#' \emph{Biological Conservation} \strong{61}: 1–10.
#' @examples
#' library(ape)
#' library(Matrix)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
#'   c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
#'   dimnames = list(paste0("g", 1:6), tree$tip.label))
#'
#' PD(com, tree)
#' @rdname PD
#' @export
PD <- function(x, phy){
  if(!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  if(length(setdiff(colnames(x), phy$tip.label)) > 0)
    stop("There are species labels in community matrix missing in the tree!")
  if(length(setdiff(phy$tip.label, colnames(x))) > 0)
    phy <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
  x <- x[,intersect(phy$tip.label, colnames(x))]
  z <- phylo_community(x, phy)
  (z$Matrix %*% z$edge.length)[,1]
}


