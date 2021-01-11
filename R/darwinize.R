#' Darwin's naturalization conundrum
#'
#' \code{darwinize} computes standard effective size of
#' mean phylogenetic distance between nearest neighbors.
#'
#' @param x A data frame or vector containing a species list
#' @param phy A phylogenetic tree
#' @param iter Numeric, the number of permutations
#' @rdname darwinize
#' @keywords bioregion
#' @importFrom ape keep.tip
#' @importFrom stats cophenetic as.dist sd
#' @return a vector with the mean distance between species
#'
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#' @examples
#' require(ape)
#' data(africa)
#' tree <- africa$phylo
#' set.seed(1809)
#' y <- sample(tree$tip.label, 10, replace = TRUE)
#'
#' darwinize(x = y, phy = tree, iter = 100)
#'
#' @export
darwinize <- function(x, phy, iter){
  x <- unique(as.character(x))
  d <- keep.tip(phy, x)
  obs <- as.dist(cophenetic(d))
  sub_tree <- phy

  y <- NULL
  for(i in seq_len(iter)){
    sub_tree$tip.label <- sub_tree$tip.label[sample(length(sub_tree$tip.label))]
    v <- keep.tip(sub_tree, x)
    v <- as.dist(cophenetic(v))
    y[i] <- mean(v)
  }
  res <- (mean(obs)-mean(y))/sd(y)
  res
}
