#' Convert community data to sparse matrix.
#'
#' This function converts community data frame to compressed sparse matrix. 
#' It is especially developed for handling big data.
#'
#' @param dat A community data frame with two columns, grids and species 
#' (or taxa)
#' @param method one of "phylo" (the default) corresponding to phylogenetic 
#' beta diversity, or "nonphylo" for beta diversity.
#' @rdname sampl2sparse
#' @keywords bioregion
#' @importFrom Matrix sparseMatrix
#' @importFrom data.table as.data.table
#'
#' @return A compressed sparse community matrix of sites by species
#'
#' @examples
#' require(data.table)
#' fdir <- system.file("NGAplants", package="phyloregion")
#' files <- file.path(fdir, dir(fdir))
#' x <- data.frame(raster2comm(files))
#' comm <- sampl2sparse(x)
#' @export sampl2sparse
sampl2sparse <- function(dat,
                         grids="grids",
                         species="species",
                         method="phylo", 
                         verbose = TRUE){
  
  if (verbose) {
    message("Converting community data to sparse matrix")
  }
  
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
