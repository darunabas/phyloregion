
#' Phyloregions for functional traits and phylogeny
#'
#' Generates a sparse community matrix as input for clustering regions
#' based on the similairity of functional traits across species.
#' @param x A community data in long format with one column representing
#' sites labeled \dQuote{grids} and another column representing species
#' labeled \dQuote{species}.
#' @param trait A data frame or matrix object with the first column
#' labeled \dQuote{species} containing the taxonomic groups to be
#' evaluated whereas the remaining columns have the various functional
#' traits. The variables must be a mix of numeric and categorical values.
#' @param cut The slice time.
#' @param phy is a dated phylogenetic tree with branch lengths stored
#' as a phylo object (as in the ape package).
#' @param bin The desired number of clusters or bins.
#' @param na.rm Logical, whether NA values should be removed or not.
#' @rdname counts
#' @keywords bioregion
#' @importFrom clustMixType kproto
#' @importFrom Matrix Matrix
#' @importFrom ape keep.tip
#' @return Function returns a community data frame that captures the
#' count of each species based on its cluster membership.
#' @export
counts <- function (x, trait, cut = NULL, phy = NULL, bin = 10, na.rm = FALSE)
{

  x$species <- gsub(" ", "_", x$species)

  if (!is.null(phy)) {
    subphy <- keep.tip(phy, intersect(phy$tip.label, x$species))
    submat <- subset(x, x$species %in% intersect(phy$tip.label, x$species))
    tx <- get_clades(subphy, cut = cut, k = bin)
    z <- length(tx)
    memb <- rep(seq_len(z), lengths(tx))
    sp <- sparseMatrix(seq_along(memb), j = memb, dims=c(length(memb),z),
                       dimnames=list(unlist(tx), as.character(seq_len(z))))
  } else {
    trait <- as.data.frame(trait)
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    index <- intersect(x$species, trait$species)
    trait <- subset(trait, trait$species %in% index)

    row.names(trait) <- trait$species
    trait <- trait[ , !colnames(trait) %in% "species"]

    x1 <- Filter(is.character, trait)
    x1[] <- lapply(x1, as.factor)

    x2 <- Filter(is.numeric, trait)
    x2[] <- apply(x2, 2, as.numeric)
    m <- data.frame(x1, x2)

    # apply k-prototypes
    m <- kproto(m, bin = bin, na.rm = na.rm)

    memb <- m$cluster
    #    names(memb) <- labels(g1)
    #    memb <- memb[!(memb==0)]
    z <- length(unique(memb))

    sp <- sparseMatrix(seq_along(memb), j = memb, dims=c(length(memb),z),
                       dimnames=list(labels(memb), as.character(seq_len(z))))

    submat <- subset(x, x$species %in% index)
  }

  M <- long2sparse(submat)

  mx <- M %*% sp
  raw_mtx <- mx/rowSums(mx)
  res <- sparse2long(mx)
  return(list(comm_dat = res, bins = z, raw_mtx = raw_mtx))
}
