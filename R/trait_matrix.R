#' Phyloregions for functional traits
#'
#' Generates a sparse community matrix as input for clustering regions
#' based on the similairity of functional traits across species.
#' @param x A community data in long format with one column representing
#' sites labeled \dQuote{grids} and another column representing species
#' labeled \dQuote{species}.
#' @param trait A data frame or matrix object with the first column
#' labeled \dQuote{species} containing the texonomic groups to be
#' evaluated whereas the remaining columns have the various functional
#' traits. The variables must be a mix of numeric and categorical values.
#' @param num A vector of column index of the numerical trait variable.
#' @param bin A vector of column index of the binary trait variable.
#' @param cat A vector of column index of the categorical trait variable.
#' @param cut The slice time.
#' @param phy is a dated phylogenetic tree with branch lengths stored
#' as a phylo object (as in the ape package).
#' @param k The desired number of clusters for single numeric trait
#' variables.
#' @param method Whether to compute phylogenetic
#' (method = \dQuote{phylo}), functional (method = \dQuote{trait}),
#' single functional categorical variables (method = \dQuote{single_cat}),
#' or for single numeric variables (method = \dQuote{single_num}).
#' @rdname trait_matrix
#' @keywords bioregion
#' @importFrom dbscan dbscan
#' @importFrom stats as.dist setNames
#' @importFrom utils stack
#' @importFrom Matrix Matrix
#' @importFrom ape keep.tip
#' @return Function returns a community data frame that captures the
#' count of each species based on its cluster membership.
#' @export
trait_matrix <- function (x, trait, num = NULL, bin = NULL,
                          cat = NULL, cut = NULL, phy = NULL,
                          method="trait", k = NULL)
{

  x$species <- gsub(" ", "_", x$species)

  if (method=="phylo") {
    subphy <- keep.tip(phy, intersect(phy$tip.label, x$species))
    submat <- subset(x, x$species %in% intersect(phy$tip.label, x$species))
    tx <- get_clades(subphy, cut = cut)
    z <- length(tx)
    memb <- rep(seq_len(z), lengths(tx))
    names(memb) <- unlist(tx)

  } else if(method=="trait") {
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    row.names(trait) <- trait$species
    g <- gower_distance(trait, method = "gower",
                        idnum=num, idbin=bin, idcat=cat)
    gc()
    g1 <- as.dist(g)
    g2 <- dbscan(g1, minPts = 5, eps = 0.05)
    gc()

    memb <- g2$cluster
    names(memb) <- labels(g1)
    memb <- memb[!(memb==0)]
    z <- length(unique(memb))
    submat <- subset(x, x$species %in% intersect(x$species, names(memb)))
  } else if (method=="single_cat"){
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    row.names(trait) <- trait$species
    zz <- trait[, cat, drop=FALSE]
    names(zz) <- "foo"

    zz$foo[zz$foo==""] <- "noise"
    memb <- as.numeric(factor(zz[,1]))
    names(memb) <- row.names(zz)
    z <- length(unique(memb))
    submat <- subset(x, x$species %in% intersect(x$species, names(memb)))
  } else if (method=="single_num"){
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    row.names(trait) <- trait$species
    zz <- trait[, num, drop=FALSE]
    names(zz) <- "foo"
    zz$cluster <- findInterval(zz$foo,
                               unique(quantile(zz$foo, probs = seq(0, 1, 1/k))),
                               rightmost.closed = TRUE)
    memb <- zz$cluster
    names(memb) <- row.names(zz)
    z <- length(unique(memb))
    submat <- subset(x, x$species %in% intersect(x$species, names(memb)))
  }

  M <- long2sparse(submat)
  mx <- Matrix(0, dim(M)[[1]], z)
  tmp <- Matrix(0, dim(M)[[1]], dim(M)[[2]])
  rownames(mx) <- rownames(M)
  rownames(tmp) <- rownames(M)
  colnames(tmp) <- names(memb)[order(memb, decreasing = FALSE)]
  colnames(mx) <- colnames(M)[1:z]

  for (i in 1:dim(M)[[1]]) {
    tmp[i, ] <- as.numeric(M[i, names(memb)[order(memb, decreasing = FALSE)]])
    for (j in 1:z) {
      names <- names(memb)[memb == j]
      mx[i, j] <- sum(tmp[i, names])
      colnames(mx)[j] <- names[[1]]
    }
  }

  res <- sparse2long(mx)
  return(list(comm_dat = res, k = z))
}
