#' Fits Grade of membership models for biogeographic regionalization
#'
#' Generates grade of membership, “admixture”, “topic”
#' or “Latent Dirichlet Allocation” models, by representing sampling units
#' as partial memberships in multiple groups. It can group regions
#' based on phylogenetic information or functional traits.
#'
#' Mapping phylogenetic regions (phyloregions) involves successively slicing
#' the phylogenetic tree at various time depths (e.g., from 1, 2, 3, 4, to 5
#' million years ago (Ma)), collapsing nodes and ranges that originated at
#' each time depth, and generating a new community matrix based on the presence
#' or absence of each lineage in a grid cell. A grade of membership model is
#' then fitted to the reduced community matrix. To map functional trait regions
#' (traitregions), the function uses k-means to cluster species based on their
#' functional traits, often for mixed-type data including categorical and
#' numeric functional traits. The ranges for each species in each resulting
#' cluster are collapsed to generate a new community matrix based on the
#' presence or absence of cluster representative in a grid cell. A grade of
#' membership model is then fitted to the new reduced community matrix. Mapping
#' bioregions for taxonomic diversity is based on fitting a grade of
#' membership model directly to the original community matrix that is often
#' represented with species in the columns and sites as rows.
#'
#' @param x A community data in long format with one column representing
#' sites labeled \dQuote{grids} and another column representing species
#' labeled \dQuote{species}.
#' @param trait A data frame or matrix object with the first column
#' labeled \dQuote{species} containing the taxonomic groups to be
#' evaluated whereas the remaining columns have the various functional
#' traits. The variables must be a mix of numeric and categorical values.
#' @param phy is a dated phylogenetic tree with branch lengths stored
#' as a phylo object (as in the ape package).
#' @param cut The slice time for the phylogenetic tree.
#' @param bin The desired number of clusters or bins.
#' @param na.rm Logical, whether NA values should be removed or not.
#' @param K The number of latent topics. If length(K)>1,
#' topics will find the Bayes factor (vs a null single topic model)
#' for each element and return parameter estimates for the highest
#' probability K.
#' @param shape Optional argument to specify the Dirichlet prior
#' concentration parameter as shape for topic-phrase probabilities.
#' Defaults to 1/(K*ncol(counts)). For fixed single K, this can also
#' be a ncol(counts) by K matrix of unique shapes for each topic element.
#' @param initopics Optional start-location for
#' \eqn{[\theta_1, \ldots, \theta_K]}, the
#' topic-phrase probabilities. Dimensions must accord with the smallest
#' element of K. If NULL, the initial estimates are built by incrementally
#' adding topics.
#' @param tol An indicator for whether or not to calculate the Bayes factor
#' for univariate K. If length(K)>1, this is ignored and Bayes factors are
#' always calculated.
#' @param bf An indicator for whether or not to calculate the Bayes factor
#' for univariate K. If length(K)>1, this is ignored and Bayes factors are
#' always calculated.
#' @param kill For choosing from multiple K numbers of topics (evaluated
#' in increasing order), the search will stop after kill consecutive drops
#' in the corresponding Bayes factor. Specify kill=0 if you want Bayes factors
#' for all elements of K.
#' @param ord If \code{TRUE}, the returned topics (columns of theta) will
#' be ordered by decreasing usage (i.e., by decreasing \code{colSums(omega)}).
#' @param verb A switch for controlling printed output. verb > 0 will
#' print something, with the level of detail increasing with verb.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname fitgom
#' @keywords bioregion
#' @importFrom clustMixType kproto
#' @importFrom Matrix Matrix
#' @importFrom ape keep.tip
#' @importFrom maptpx topics
#' @return An topics object list with entries
#' \itemize{
#'   \item \code{K} The number of latent topics estimated. If input
#'   \code{length(K)>1}, on output this is a single value corresponding to
#'   the model with the highest Bayes factor.
#'   \item \code{theta} The ncol{counts} by K matrix of estimated
#'   topic-phrase probabilities.
#'   \item \code{omega} The nrow{counts} by K matrix of estimated
#'   document-topic weights.
#'   \item \code{BF} The log Bayes factor for each number of topics
#'   in the input K, against a null single topic model.
#'   \item \code{D} Residual dispersion: for each element of K, estimated
#'   dispersion parameter (which should be near one for the multinomial),
#'   degrees of freedom, and p-value for a test of whether the true dispersion
#'   is >1.
#'   \item \code{X} The input community matrix as a sparse matrix.
#' }
#' @examples
#' library(terra)
#' data(africa)
#' names(africa)
#' p <- vect(system.file("ex/sa.json", package = "phyloregion"))
#' m <- fitgom(x=sparse2long(africa$comm), K=3)
#'
#' COLRS <- phyloregion:::hue(m$K)
#' plot_pie(m$omega, pol = p, col=COLRS)
#' @export
fitgom <- function (x, trait = NULL, cut = NULL, phy = NULL, bin = 10,
                    na.rm = FALSE, K, shape = NULL, initopics = NULL, tol = 0.1,
                    bf = TRUE, kill = 2, ord = TRUE, verb = 1, ...)
{

  if(is(x, "sparseMatrix")) {
    x <- sparse2long(x)
  }
  x$species <- gsub(" ", "_", x$species)

  if (!is.null(phy)) {
    subphy <- keep.tip(phy, intersect(phy$tip.label, x$species))
    submat <- subset(x, x$species %in% intersect(phy$tip.label, x$species))
    tx <- get_clades(subphy, cut = cut, k = bin)

    # SPECIES THAT CAN BE LUMPED
    m <- do.call(rbind, lapply(tx, function(y) cbind(species=y,
                                                     cluster=parent.frame()$i[],
                                                     grouping="lumper",
                                                     ntax=length(y))))
    r <- data.frame(m[m[, "ntax"] > 1,])
    # MASTER_DISTINCT -- SPECIES WHICH ARE DISTINCT AND CAN'T BE LUMPED
    w <- submat[!submat$species %in% as.vector(r$species), ]
    spp <- unique(as.character(r$cluster))

    spo <- lapply(spp, function(y) {
      v <- r[r$cluster==y, ]
      v <- submat[submat$species %in% as.vector(v$species), ]
      cbind(v["grids"], species=v$species[1])
    })

    z <- rbind(w, do.call(rbind, spo))
    mx <- long2sparse(z)
  } else if (!is.null(trait)) {
    trait <- as.data.frame(trait)
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    index <- intersect(x$species, trait$species)
    submat <- x[x$species %in% index, ]
    trait <- subset(trait, trait$species %in% index)

    row.names(trait) <- trait$species
    trait <- trait[ , !colnames(trait) %in% "species"]

    x1 <- Filter(is.character, trait)
    x1[] <- lapply(x1, as.factor)

    x2 <- Filter(is.numeric, trait)
    x2[] <- apply(x2, 2, as.numeric)
    mm <- data.frame(x1, x2)

    # apply k-prototypes
    mm <- kproto(mm, k = bin, na.rm = na.rm)
    memb <- mm$cluster

    foox <- function(x) {
      k <- max(x)
      res <- vector("list", k)
      for (i in seq_len(k)) {
        res[[i]] <- names(x)[x==i]
      }
      res
    }
    d <- foox(memb)
    m <- do.call(rbind, lapply(d, function(y) cbind(species=y,
                                                    cluster=parent.frame()$i[],
                                                    grouping="lumper",
                                                    ntax=length(y))))
    r <- data.frame(m[m[, "ntax"] > 1,])
    # MASTER_DISTINCT -- SPECIES WHICH ARE DISTINCT AND CAN'T BE LUMPED
    w <- submat[!submat$species %in% as.vector(r$species), ]
    spp <- unique(as.character(r$cluster))

    spo <- lapply(spp, function(y) {
      v <- r[r$cluster==y, ]
      v <- submat[submat$species %in% as.vector(v$species), ]
      cbind(v["grids"], species=v$species[1])
    })
    z <- rbind(w, do.call(rbind, spo))
    mx <- long2sparse(z)
  } else {
    mx <- long2sparse(x)
  }
  res <- topics(mx, K, shape = shape, initopics = initopics, tol = tol,
                bf = bf, kill = kill, ord = ord, verb = verb, ...)
  return(res)

}
