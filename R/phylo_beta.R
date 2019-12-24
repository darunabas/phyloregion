phylo_com <- function(tip, phy) {
  if (!inherits(phy, "phylo"))
    stop("object \"phy\" is not of class \"phylo\"")
  done_v <- logical(length(phy$tip.label) + phy$Nnode)
  rootnd <- getRoot(phy)
  if (is.character(tip)) tip <- match(tip, c(phy$tip.label, phy$node.label))
  tip <- as.integer(tip)
  res <- pvec <- integer(max(phy$edge))
  pvec[phy$edge[, 2]] <- phy$edge[, 1]

  res[seq_along(tip)] <- tip
  l <- length(tip) + 1L
  res[l] <- rootnd
  done_v[rootnd] <- TRUE
  for (k in tip) {
    nd <- pvec[k]
    done <- done_v[nd]
    while (!done) {
      done_v[nd] <- TRUE
      l <- l + 1
      res[l] <- nd
      nd <- pvec[nd]
      done <- done_v[nd]
    }
  }
  sort(res[1:l])
}


phylo_community <- function(x, phy) {
  el <- numeric(max(phy$edge))
  el[phy$edge[, 2]] <- phy$edge.length
  x <- x[, phy$tip.label]
  if (is.matrix(x) | is(x, "sparseMatrix")) {
    x <- as.splits(x)
  }
  if (is.character(x) | is.numeric(x)) y <- list(phylo_com(x, phy))
  if (is.list(x)) {
    y <- lapply(x, function(x, phy) phylo_com(x, phy), phy)
  }
  if (is.null(y)) return(NULL)
  attr(y, "edge.length") <- el
  attr(y, "labels") <- c(phy$tip.label, as.character(Ntip(phy)
  + (1:Nnode(phy))))
  class(y) <- c("phylo_community", "splits")
  y
}


#' Phylogenetic beta diversity
#'
#' \code{beta_diss} and \code{phylobeta_core} computes efficiently for
#' large community matrices and trees the necessary quantities used by the
#' betapart package to compute pairwise and multiple-site phylogenetic
#' dissimilarities.
#'
#' @aliases phylobeta phylobeta_core
#' @param x an object of class Matrix or matrix
#' @param phy a phylogenetic tree (object of class phylo)
#' @param index.family family of dissimilarity indices, partial match of
#' "sorensen" or "jaccard".
#' @keywords phyloregion
#' @seealso \code{\link{read.community}}, \code{\link{phylo.betapart.core}}
#' @examples
#' library(ape)
#' tree <- read.tree(text = "((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1, 0, 1, 1, 0, 0,
#'   1, 0, 0, 1, 1, 0,
#'   1, 1, 1, 1, 1, 1,
#'   0, 0, 1, 1, 0, 1), 6, 4,
#' dimnames = list(paste0("g", 1:6), tree$tip.label))
#' pbc <- phylobeta_core(com, tree)
#' pb <- phylobeta(com, tree)
#' @rdname phylobeta
#' @importFrom fastmatch fmatch
#' @importFrom betapart phylo.beta.multi phylo.beta.pair
#' @importFrom phangorn getRoot
#' @export
phylobeta_core <- function(x, phy) {
  if (!identical(sort(colnames(x)), sort(phy$tip.label)))
    stop("Labels of community matrix and tree differ!")
  x <- phylo_community(x, phy)
  l <- length(x)
  el <- attr(x, "edge.length")
  pd_tmp <- vapply(x, function(x, el) sum(el[x]), 0, el)

  Labels <- names(x)
  class(x) <- NULL
  SHARED <- vector("numeric", l * (l - 1) / 2)
  B <- vector("numeric", l * (l - 1) / 2)
  C <- vector("numeric", l * (l - 1) / 2)

  k <- 1
  for (i in 1:(l - 1)) {
    xi <- x[[i]]
    for (j in (i + 1):l) {
      #   sum(el[fast_intersect(x[[i]], x[[j]])])
      SHARED[k] <- sum(el[xi[fmatch(x[[j]], xi, 0L)]])
      B[k] <- pd_tmp[i] - SHARED[k]
      C[k] <- pd_tmp[j] - SHARED[k]
      k <- k + 1
    }
  }

  sum.not.shared <- B + C
  max.not.shared <- pmax(B, C)
  min.not.shared <- pmin(B, C)

  at <- structure(list(Labels = Labels, Size = l, class = "dist", Diag = FALSE,
    Upper = FALSE), .Names = c("Labels", "Size", "class", "Diag", "Upper"))
  attributes(SHARED) <- at
  attributes(sum.not.shared) <- at
  attributes(max.not.shared) <- at
  attributes(min.not.shared) <- at
  res <- list(sumSi = sum(pd_tmp), St = sum(el), shared = SHARED,
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
  if (index.family == "simpson") {
    z <- p[[1]]
  }
  else if (index.family == "sorensen") {
    z <- p[[1]]
  }
  else if (index.family == "jaccard") {
    z <- p[[3]]
  }
  return(z)
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
#' @keywords bioregion
#' @examples
#' data(africa)
#' tree <- africa$phylo
#' x <- sampl2sparse(africa$comm)
#'
#' subphy <- match_phylo_comm(tree, x)$phy
#' submat <- match_phylo_comm(tree, x)$com
#' @importFrom ape keep.tip
#' @export
match_phylo_comm <- function(phy, comm) {
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
  res$comm <- comm[, index]
  res$phy <- keep.tip(phy, index)
  return(res)
}

#' Taxonomic (non-phylogenetic) beta diversity
#' @rdname beta_diss
#' @param x an object of class Matrix or matrix.
#' @param index.family family of dissimilarity indices, partial match of
#' "sorensen" or "jaccard".
#' @importFrom Matrix Matrix tcrossprod colSums
#' @importFrom betapart beta.pair beta.multi
#' @seealso \code{\link{betapart.core}}
#' @export
## non phylogenetic version
beta_core <- function(x) {
  if (!inherits(x, "Matrix")) x <- Matrix(x)
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
  p <- beta.pair(computations, index.family = index.family)
  if (index.family == "simpson") {
    res <- p[[1]]
  }
  else if (index.family == "sorensen") {
    res <- p[[3]]
  }
  else if (index.family == "jaccard") {
    res <- p[[3]]
  }
  return(res)
}
