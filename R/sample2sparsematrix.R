#' Convert community data to sparse matrix.
#'
#' This function converts community data frame to compressed sparse matrix.
#' It is fast and especially developed for handling large datasets.
#'
#' @param dat A community data frame with at least two columns, grids and species
#' @param grids column name of the column containing grid cells
#' @param species column name of the column containing  the species / taxa names
#' @param method one of "phylo" (the default) corresponding to phylogenetic
#' beta diversity, or "nonphylo" for beta diversity.
#' @rdname sampl2sparse
#' @importFrom Matrix sparseMatrix
#' @importFrom data.table as.data.table
#'
#' @return A compressed sparse community matrix of sites by species
#'
#' @examples
#' data(africa)
#' head(africa$comm)
#' comm <- sampl2sparse(africa$comm)
#' @export sampl2sparse
sampl2sparse <- function(dat, grids="grids", species="species", method="phylo"){
  dat <- dat[, c(grids, species)]
  names(dat) <- c("grids", "species")
  grids <- factor(as.character(dat$grids))
  species <- factor(as.character(dat$species))
  if (method == "phylo") {
    res <- sparseMatrix(as.integer(grids), as.integer(species),
                        dimnames = list(levels(grids), levels(species)))
  }
  else if (method == "nonphylo") {
    res <- sparseMatrix(as.integer(grids), as.integer(species),
                      x=rep(1L, length(grids)),
                      dimnames = list(levels(grids), levels(species)))
  }
  res
}
