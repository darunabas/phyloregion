#' Subset trees from posterior distribution of trees.
#'
#' This function randomly samples a subset of trees from a posterior 
#' distribution of trees derived from multiple runs of MrBayes.
#'
#' @param wd A path to the working directory with the distributions of 
#' multiple phylogenetic trees.
#' @param n The desired number of subsets of trees. This defaults to 100.
#' @param \dots arguments passed among methods.
#' 
#' @rdname tree.sampler
#' @importFrom ape read.tree
#' @export
#'
#' @return An object of class "multiPhylo" with the subset of trees.
#'
#' @author Dominic Bennett & Harith Farooq \email{harithmorgadinho@gmail.com}
#'
#'
#' @export


tree.sampler <- function(wd, n = 100, ...) {
  read_treelines <- function(fl) {
    readLines(con = fl)
  }
  sampler <- function(i) {
    read.tree(text = sample(tree_lines, 1))
  }
  fls <- list.files(path = wd, pattern = '.trees')
  fls <- file.path(wd, fls)
  tree_lines <- lapply(X = fls, FUN = read_treelines)
  tree_lines <- unlist(tree_lines)
  trees <- lapply(X = seq_len(n), FUN = sampler)
  class(trees) <- 'multiPhylo'
  trees
}

