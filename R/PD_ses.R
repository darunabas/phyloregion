rt <- function (phy) {
    phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
    return(phy)
}

#' Phylogenetic diversity standardized for species richness
#'
#' This function computes the standard effect size of PD by correcting for
#' changes in species richness. The novelty of this function is its ability
#' to utilize sparse community matrix making it possible to
#' efficiently randomize very large community matrices spanning thousands of
#' taxa and sites.
#'
#' @param x  a (sparse) community matrix, i.e. an object of class matrix or Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @param reps Number of replications.
#' @param model The null model for separating patterns from processes and
#' for contrasting against alternative hypotheses. Available null models
#' include:
#' \itemize{
#' \item \dQuote{tipshuffle}: shuffles tip labels multiple times.
#' \item \dQuote{rowwise}: shuffles sites (i.e., varying richness) and
#' keeping species occurrence frequency constant.
#' \item \dQuote{colwise}: shuffles species occurrence frequency and
#' keeping site richness constant.}
#'
#' @param \dots Further arguments passed to or from other methods.
#' @rdname PD_ses
#' @importFrom stats sd var
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
#' PD_ses(com, tree, model="rowwise")
#'
#' @return A data frame of results for each community or grid cell
#' \itemize{
#'   \item grids: Site identity
#'   \item richness: Number of taxa in community
#'   \item PD_obs: Observed PD in community
#'   \item pd_rand_mean: Mean PD in null communities
#'   \item pd_rand_sd: Standard deviation of PD in null communities
#'   \item pd_obs_rank: Rank of observed PD vs. null communities
#'   \item pd_obs_z: Standardized effect size of PD vs. null communities \eqn{= (PD_obs - pd_rand_mean) / pd_rand_sd}
#'   \item pd_obs_p: P-value (quantile) of observed PD vs. null communities \eqn{= mpd_obs_rank / iter + 1}
#'   \item reps: Number of replicates
#' }
#' @export
PD_ses <- function(x, phy,
                   model = c("tipshuffle", "rowwise", "colwise"),
                   reps = 1000, ...) {

    colnames(x) <- gsub(" ", "_", colnames(x))
    p <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
    x <- x[, intersect(p$tip.label, colnames(x))]

    PD_obs <- PD(x, p)
    pd.rand <- switch(model,
                      tipshuffle = lapply(seq_len(reps), function(i)
                          PD(x, rt(p))),
                      rowwise = lapply(seq_len(reps), function(i)
                          PD(x[sample(nrow(x)),], p)),
                      colwise = lapply(seq_len(reps), function(i)
                          PD(x[,sample(ncol(x))], p)))

    y <- do.call(rbind, pd.rand)
    pd_rand_mean <- apply(X = y, MARGIN = 2, FUN = mean, na.rm = TRUE)
    pd_rand_sd <- apply(X = y, MARGIN = 2, FUN = var, na.rm = TRUE)
    #pd_rand_sd[pd_rand_sd == 0] <- 1
    zscore <- (PD_obs - pd_rand_mean)/sqrt(pd_rand_sd)
    pd_obs_rank <- apply(X = rbind(PD_obs, y), MARGIN = 2, FUN = rank)[1, ]
    pd_obs_rank <- ifelse(is.na(pd_rand_mean), NA, pd_obs_rank)

    m <- data.frame(grids=rownames(x), PD_obs, pd_rand_mean, pd_rand_sd,
                    pd_obs_rank, zscore,
                    pd_obs_p = pd_obs_rank/(reps + 1), reps = reps,
                    row.names = row.names(x))


    z <- data.frame(table(sparse2long(x)$grids))
    names(z) <- c("grids", "richness")
    res <- Reduce(function(x, y) merge(x, y, by="grids",all=TRUE) , list(z, m))
    res
}
