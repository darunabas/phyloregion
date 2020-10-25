rt <- function (phy) {
    phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
    return(phy)
}

#' Phylogenetic beta diversity standardized for species beta diversity
#'
#' This function computes the standard effect size of phylogenetic beta
#' diversity by correcting for changes in species beta diversity. The
#' novelty of this function is its ability to utilize sparse community
#' matrix making it possible to efficiently randomize very large
#' community matrices spanning thousands of taxa and sites.
#'
#' @param x  a (sparse) community matrix, i.e., an object of class
#' matrix or Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @param reps Number of replications.
#' @param index.family the family of dissimilarity indices including
#' \dQuote{simpson}, \dQuote{sorensen} and \dQuote{jaccard}.
#' @param model The null model for separating patterns from processes and
#' for contrasting against alternative hypotheses. Available null models
#' include:
#' \itemize{
#' \item \dQuote{tipshuffle}: shuffles phylogenetic tip labels multiple times.
#' \item \dQuote{rowwise}: shuffles sites (i.e., varying richness) and
#' keeping species occurrence frequency constant.
#' \item \dQuote{colwise}: shuffles species occurrence frequency and
#' keeping site richness constant.}
#'
#' @param \dots Further arguments passed to or from other methods.
#' @rdname phylobeta_ses
#' @importFrom stats sd
#' @importFrom ape keep.tip
#' @references
#' Proches, S., Wilson, J.R.U. & Cowling, R.M. (2006) How much evolutionary
#' history in a 10 x 10m plot? \emph{Proceedings of Royal Society B}
#' \strong{273}: 1143-1148.
#' @examples
#' library(ape)
#' library(Matrix)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
#'   c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
#'   dimnames = list(paste0("g", 1:6), tree$tip.label))
#'
#' phylobeta_ses(com, tree, model="rowwise")
#'
#' @return A data frame of results for each community or grid cell
#' \itemize{
#'   \item phylobeta_obs: Observed phylobeta in community
#'   \item phylobeta_rand_mean: Mean phylobeta in null communities
#'   \item phylobeta_rand_sd: Standard deviation of phylobeta in null communities
#'   \item phylobeta_obs_z: Standardized effect size of phylobeta vs.
#'   null communities \eqn{= (phylobeta_obs - phylobeta_rand_mean)
#'   / phylobeta_rand_sd}
#'   \item reps: Number of replicates
#' }
#' @export
phylobeta_ses <- function(x, phy, index.family="simpson",
                   model = c("tipshuffle", "rowwise", "colwise"),
                   reps = 1000, ...) {

    colnames(x) <- gsub(" ", "_", colnames(x))
    p <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
    x <- x[, intersect(p$tip.label, colnames(x))]

    new_phylobeta <- function(mat, tree, index.family) {
        switch(index.family,
               simpson = phylobeta(mat, tree, index.family = "sorensen")[[1]],
               sorensen = phylobeta(mat, tree, index.family= "sorensen")[[3]],
               jaccard = phylobeta(mat, tree, index.family = "jaccard")[[3]])
    }

    pbd_obs <- new_phylobeta(x, p, index.family)
    mean_obs_z <- pbd_obs


    dnam <- labels(pbd_obs)
    sort_nam <- function(x, dnam) {
        d <- as.matrix(x)[dnam,dnam]
        as.dist(d)
    }

    pbd.rand <- switch(model,
                       tipshuffle = sapply(seq_len(reps), function(i)
                           as.vector(new_phylobeta(x, rt(p), index.family))),
                       rowwise = sapply(seq_len(reps), function(i)
                           as.vector(sort_nam(new_phylobeta(x[sample(nrow(x)),],
                                                        p, index.family), dnam))),
                       colwise = sapply(seq_len(reps), function(i)
                           as.vector(new_phylobeta(x[,sample(ncol(x))],
                                                   p, index.family))))

    pbd_rand_mean <- apply(X = pbd.rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
    mean_obs_z[] <- pbd_rand_mean
    pbd_rand_mean <- mean_obs_z

    pbd_rand_sd <- apply(X = pbd.rand, MARGIN = 1, FUN = sd, na.rm = TRUE)
    pbd_rand_sd[pbd_rand_sd == 0] <- 1
    mean_obs_z[] <- pbd_rand_sd
    pbd_rand_sd <- mean_obs_z

    pbd_obs_z <- (as.vector(pbd_obs) - pbd_rand_mean)/pbd_rand_sd
    mean_obs_z[] <- pbd_obs_z
    pbd_obs_z <- mean_obs_z

    return(list(phylobeta_obs = pbd_obs,
                phylobeta_rand_mean = pbd_rand_mean,
                phylobeta_rand_sd = pbd_rand_sd,
                phylobeta_obs_z = pbd_obs_z, reps = reps))

}
