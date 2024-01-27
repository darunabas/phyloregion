#' Phylogenetic diversity
#'
#' \code{PD} calculates Faith's (1992) phylogenetic diversity and relative
#' phylogenetic diversity.
#'
#' @param x a community matrix, i.e. an object of class matrix or Matrix or an
#' object of class phyloseq.
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
#' # Relative phylogenetic diversity
#' RPD(com, tree)
#' @rdname PD
#' @export
PD <- function(x, phy){
  if(inherits(x, "phyloseq")){
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      if(missing(phy)) phy <- phyloseq::phy_tree(x)
      otu <- as(phyloseq::otu_table(x), "matrix")
      if (phyloseq::taxa_are_rows(x)) otu <- t(otu)
      x <- Matrix(otu, sparse = TRUE)
    }
  }
  if(inherits(x, "matrix") && ncol(x)>2) x <- Matrix(x, sparse = TRUE)
  if(!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  if(length(setdiff(colnames(x), phy$tip.label)) > 0)
    stop("There are species labels in community matrix missing in the tree!")
  if(length(setdiff(phy$tip.label, colnames(x))) > 0)
    phy <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
  x <- x[,intersect(phy$tip.label, colnames(x))]
  z <- phylo_community(x, phy)
  (z$Matrix %*% z$edge.length)[,1]
}

#' @rdname PD
#' @export
RPD <- function(x, phy) {
  if(inherits(x, "matrix") && ncol(x)>2) x <- Matrix(x, sparse = TRUE)
  if(!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  if(length(setdiff(colnames(x), phy$tip.label)) > 0)
    stop("There are species labels in community matrix missing in the tree!")
  dat <- match_phylo_comm(phy, x)
  phy <- dat$phy
  phy_alt <- phy
  # convert **non-zero** branch lengths to same value (1)
  nb <- unlist(lapply(phy_alt$edge.length, function(x) !identical(x, 0)))
  phy_alt$edge.length[nb] <- rep(x=1, times=length(phy_alt$edge.length[nb])) 
  phy_alt$edge.length <- phy_alt$edge.length / sum(phy_alt$edge.length)
  # rescale original phy so total length is 1
  phy$edge.length <- phy$edge.length / sum(phy$edge.length)
  # 1. PD
  pd <- PD(x, phy)
  pd <- data.frame(grids=names(pd), pd)
  # 2. PD_ALT
  pd_alt <- PD(x, phy_alt)
  pd_alt <- data.frame(grids=names(pd_alt), pd_alt)
  w <- Reduce(function(x, y) merge(x, y, by="grids"), list(pd, pd_alt))
  w$rpd <- w$pd / w$pd_alt
  rpd <- w$rpd
  names(rpd) <- w$grids
  return(rpd)
}

