rand_tip <- function (phy) {
    phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
    return(phy)
}

#' PD standardized for species richness
#'
#' This function computes the standard effect size of PD by correcting for
#' changes in species richness
#'
#' @param x  a community matrix, i.e. an object of class matrix or Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @param iter Number of iterations.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname PD_ses
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
#' PD_ses(com, tree)
#'
#' @return A data frame of results for each community or grid cell
#' \itemize{
#'   \item grids: Site identity
#'   \item ntaxa: Number of taxa in community
#'   \item PD_obs: Observed PD in community
#'   \item pd_rand_mean: Mean PD in null communities
#'   \item pd_rand_sd: Standard deviation of PD in null communities
#'   \item pd_obs_rank: Rank of observed PD vs. null communities
#'   \item pd_obs_z: Standardized effect size of PD vs. null communities \eqn{= (PD_obs - pd_rand_mean) / pd_rand_sd}
#'   \item pd_obs_p: P-value (quantile) of observed PD vs. null communities \eqn{= mpd_obs_rank / iter + 1}
#'   \item iter: Number of iterations
#' }
#' @export
PD_ses <- function(x, phy, iter = 1000, ...){

    phy <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
    x <- x[, intersect(phy$tip.label, colnames(x))]

    PD_obs <- as.vector(data.frame(PD(x, phy))$PD)
    pd_rand <- t(replicate(iter, as.vector(data.frame(PD(x, rand_tip(phy)))$PD)))
    pd_rand_mean <- apply(X = pd_rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    pd_rand_sd <- apply(X = pd_rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    pd_obs_z <- (PD_obs - pd_rand_mean)/pd_rand_sd
    pd_obs_rank <- apply(X = rbind(PD_obs, pd_rand), MARGIN = 2, FUN = rank)[1, ]
    pd_obs_rank <- ifelse(is.na(pd_rand_mean), NA, pd_obs_rank)

    m <- data.frame(PD_obs, pd_rand_mean, pd_rand_sd, pd_obs_rank, pd_obs_z,
                     pd_obs_p = pd_obs_rank/(iter + 1), iter = iter,
                     row.names = row.names(x))
    m <- data.frame(grids = row.names(m), m)

    z <- data.frame(table(sparse2long(x)$grids))
    names(z) <- c("grids", "ntaxa")
    res <- Reduce(function(x, y) merge(x, y, by="grids",all=TRUE) , list(z, m))
    res
}
