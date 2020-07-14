# internally used in pd, phylobeta_core and phylo_endemism
phylo_community <- function(x, phy) {
  el <- numeric(max(phy$edge))
  el[phy$edge[, 2]] <- phy$edge.length
  x <- x[, phy$tip.label]
  anc <- Ancestors(phy, seq_along(phy$tip.label))
  anc <-  mapply(c, seq_along(phy$tip.label), anc, SIMPLIFY=FALSE)
  M <- Matrix::sparseMatrix(as.integer(rep(seq_along(anc), lengths(anc))),
                            as.integer(unlist(anc)), x = 1L)
  commphylo <- x %*% M
  commphylo@x[commphylo@x > 1e-8] <- 1
  list(Matrix = commphylo, edge.length = el)
}


#' Phylogenetic beta diversity
#'
#' \code{phylobeta_core} computes efficiently for large community matrices and
#' trees the necessary quantities used by the betapart package to compute
#' pairwise and multiple-site phylogenetic dissimilarities.
#'
#' @aliases phylobeta phylobeta_core
#' @param x an object of class Matrix or matrix
#' @param phy a phylogenetic tree (object of class phylo)
#' @param index.family family of dissimilarity indices, partial match of
#' "sorensen" or "jaccard".
#' @return \code{phylobeta_core} returns an object of class "phylo.betapart",
#' see \code{\link{phylo.betapart.core}} for details. This object can be called
#' by \code{\link{phylo.beta.pair}} or \code{\link{phylo.beta.multi}}.
#'
#' \code{phylobeta} returns a list with three phylogenetic dissimilarity
#' matrices. See \code{\link{phylo.beta.pair}} for details.
#' @keywords phyloregion
#' @seealso \code{\link{read.community}}, \code{\link{phylo.betapart.core}},
#' \code{\link{beta_core}}
#' @examples
#' library(ape)
#' tree <- read.tree(text = "((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
#'   c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
#'   dimnames = list(paste0("g", 1:6), tree$tip.label))
#' com
#'
#' pbc <- phylobeta_core(com, tree)
#' pb <- phylobeta(com, tree)
#' @rdname phylobeta
#' @author Klaus Schliep
#' @importFrom betapart phylo.beta.multi phylo.beta.pair
#' @importFrom phangorn getRoot
#' @importFrom Matrix tril tcrossprod Diagonal
#' @export
phylobeta_core <- function(x, phy) {
  if(!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  if (!identical(sort(colnames(x)), sort(phy$tip.label)))
    stop("Labels of community matrix and tree differ!")
  Labels <- rownames(x)
  x <- phylo_community(x, phy)
  pd_tmp <- (x$Matrix %*% x$edge.length)[,1]
  l <- length(pd_tmp)
  m <- l - 1L

  SHARED <- tcrossprod(x$Matrix, x$Matrix %*% Diagonal(x=x$edge.length) )
  SHARED <- as.dist(as.matrix(SHARED))
#  SHARED <- tril(SHARED, k=-1)@x

  B <- pd_tmp[rep(1:m, m:1)] - SHARED
  C <- pd_tmp[sequence(m:1) + rep(1:m, m:1)] - SHARED

  sum.not.shared <- B + C
  max.not.shared <- pmax(B, C)
  min.not.shared <- pmin(B, C)

  at <- structure(list(Labels = Labels, Size = l, class = "dist", Diag = FALSE,
        Upper = FALSE), .Names = c("Labels", "Size", "class", "Diag", "Upper"))
  attributes(SHARED) <- at
  attributes(sum.not.shared) <- at
  attributes(max.not.shared) <- at
  attributes(min.not.shared) <- at
  res <- list(sumSi = sum(pd_tmp), St = sum(x$edge.length), shared = SHARED,
              sum.not.shared = sum.not.shared,
              max.not.shared = max.not.shared, min.not.shared = min.not.shared)
  class(res) <- "phylo.betapart"
  res
}



#' @rdname phylobeta
#' @export
phylobeta <- function(x, phy, index.family = "sorensen") {
  res <- phylobeta_core(x, phy)
  p <- phylo.beta.pair(res, index.family = index.family)
  return(p)
}



#' Match taxa and in phylogeny and community matrix
#'
#' match_phylo_comm compares taxa (species, labels, tips) present in a phylogeny
#' with a community matrix. Pruning, sorting and trying to add missing species
#' on genus level if possible to match in subsequent analysis.
#'
#' Based on the function of the same name in picante but allows sparse matrices
#' and with taxa addition.
#'
#' @param phy A phylogeny
#' @param comm A (sparse) community data matrix
#' @param delete_empty_rows delete rows with no observation
#' @return A list containing the following elements, pruned and sorted to match one another:
#' \item{phy}{A phylogeny object of class phylo}
#' \item{comm}{A (sparse) community data matrix}
#' @keywords bioregion
#' @examples
#' data(africa)
#' tree <- africa$phylo
#' x <- africa$comm
#'
#' subphy <- match_phylo_comm(tree, x)$phy
#' submat <- match_phylo_comm(tree, x)$com
#' @importFrom ape keep.tip
#' @export
match_phylo_comm <- function(phy, comm, delete_empty_rows=TRUE) {
  if (!(is.data.frame(comm) | is.matrix(comm) | inherits(comm, "Matrix"))) {
    stop("Community data should be a data.frame or matrix with samples in rows
         and taxa in columns")
  }
  res <- vector("list", 2)
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to
         match phylogeny and community data")
  }
  phytaxa <- phy$tip.label
  index <- intersect(commtaxa, phytaxa)
  comm <- comm[, index, drop = FALSE]
  if(delete_empty_rows){
    comm <- comm[rowSums(comm)>0, , drop = FALSE]
  }
  phy <- keep.tip(phy, index)
  list(comm=comm, phy=phy)
}

#' Taxonomic (non-phylogenetic) beta diversity
#'
#' Data are assumed to be presence / absence (0 / 1) and all values greater zero
#' are assumed to reflect presence.
#'
#' \code{beta_core} is helper function to compute the basic quantities needed
#' for computing the "sorensen" or "jaccard" index.
#' @rdname beta_diss
#' @param x an object of class Matrix, where rows are sites and columns are
#' species.
#' @param index.family family of dissimilarity indices, partial match of
#' "sorensen" or "jaccard".
#' @return \code{beta_core} returns an object of class \code{beta_diss} like the
#' \code{\link{betapart.core}} function. This object can be called by
#' \code{\link{beta.pair}} or \code{\link{beta.multi}}.
#'
#' \code{beta_diss} returns a list with three dissimilarity matrices. See
#' \code{\link{beta.pair}} for details.
#' @importFrom Matrix Matrix tcrossprod colSums
#' @importFrom betapart beta.pair beta.multi
#' @seealso \code{\link{betapart.core}}, \code{\link{betapart}},
#' \code{\link{phylobeta}}
#' @author Klaus Schliep
#' @examples
#' data(africa)
#' x <- africa$comm
#' bc <- beta_core(x)
#' beta_sorensen <- beta_diss(x)
#' @export
## non phylogenetic version
beta_core <- function(x) {
  if (!inherits(x, "Matrix")) x <- Matrix(x, sparse=TRUE)
  x@x[x@x > 1e-8] <- 1
  shared <- as.matrix(tcrossprod(x)) # %*% t(x)
  not.shared <- abs(sweep(shared, 2, diag(shared)))
  sumSi <- sum(diag(shared))
  St <- sum(colSums(x) > 0)
  a <- sumSi - St
  sum.not.shared <- not.shared + t(not.shared)
  max.not.shared <- pmax(not.shared, t(not.shared))
  min.not.shared <- pmin(not.shared, t(not.shared))
  computations <- list(data = x, sumSi = sumSi, St = St, a = a, shared = shared,
    not.shared = not.shared, sum.not.shared = sum.not.shared,
    max.not.shared = max.not.shared,
    min.not.shared = min.not.shared)
  class(computations) <- "betapart"
  computations
}

#' @rdname beta_diss
#' @export
beta_diss <- function(x, index.family = "sorensen") {
  computations <- beta_core(x)
  res <- beta.pair(computations, index.family = index.family)
  return(res)
}
