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
#' @param x  a (sparse) community matrix, i.e. an object of class matrix or
#' Matrix.
#' @param phy a phylogenetic tree (object of class phylo).
#' @param reps Number of replications.
#' @param metric The phylodiversity measure to compute.
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
#'   \item pd_obs: Observed PD in community
#'   \item pd_rand.mean: Mean PD in null communities
#'   \item pd_rand.sd: Standard deviation of PD in null communities
#'   \item pd_obs.rank: Rank of observed PD vs. null communities
#'   \item pd_obs.z: Standardized effect size of PD vs. null communities
#'   \eqn{= (pd_obs - pd_rand.mean) / pd_rand_sd}
#'   \item pvalue: P-value (quantile) of observed PD vs. null communities
#'   \eqn{= mpd_obs_rank / iter + 1}
#'   \item reps: Number of replicates
#'   \item p_obs_c_lower: Number of times observed value < random value
#'   \item p_obs_c_upper: Number of times observed value > random value
#'   \item p_obs_p_lower: Percentage of times observed value < random value
#'   \item p_obs_p_upper: Percentage of times observed value > random value
#'   \item p_obs_q: Number of the non-NA random values used for comparison
#' }
#' @export
PD_ses <- function(x, phy,
                   model = c("tipshuffle", "rowwise", "colwise"),
                   reps = 10, metric = "pd", ...) {

    colnames(x) <- gsub(" ", "_", colnames(x))
    p <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
    x <- x[, intersect(p$tip.label, colnames(x))]

    obs <- PD(x, p)
    rand <- switch(model,
                   tipshuffle = lapply(seq_len(reps), function(i)
                       PD(x, rt(p))),
                   rowwise = lapply(seq_len(reps), function(i)
                       PD(x[sample(nrow(x)),], p)),
                   colwise = lapply(seq_len(reps), function(i)
                       PD(x[,sample(ncol(x))], p)))

    y <- do.call(rbind, rand)
    rand.mean <- apply(X = y, MARGIN = 2, FUN = mean, na.rm = TRUE)
    rand.sd <- apply(X = y, MARGIN = 2, FUN = sd, na.rm = TRUE)
    obs.z <- (obs - rand.mean)/sqrt(rand.sd)
    obs.rank <- apply(X = rbind(obs, y), MARGIN = 2, FUN = rank)[1, ]
    obs.rank <- ifelse(is.na(rand.mean), NA, obs.rank)
    pvalue <- obs.rank/(length(rand)+1)

    obs_c_upper <- numeric(length = ncol(y))
    for (i in seq_len(ncol(y))) {obs_c_upper[i] <- sum(obs[i] > y[, i])}

    obs_c_lower <- numeric(length = ncol(y))
    for (i in seq_len(ncol(y))) {obs_c_lower[i] <- sum(obs[i] < y[, i])}
    # Count the number of non-NA random values used for comparison
    obs_q <- apply(y, 2, function(x) sum(!is.na(x)))
    # Calculate p-value for upper tail
    obs_pvalue_upper <- obs_c_upper / obs_q
    # Calculate p-value for lower tail
    obs_pvalue_lower <- obs_c_lower / obs_q
    res <- data.frame(
        ID=rownames(x),
        obs, rand.mean, rand.sd, obs.z, obs.rank, obs_c_upper, obs_c_lower,
        obs_q, obs_pvalue_upper, obs_pvalue_lower, pvalue
    )

    colnames(res)[-1] <- paste(metric, colnames(res)[-1], sep = "_")
    # Calculate the significance column 'signif'
    res$signif <- with(res, ifelse(
        is.na(get(paste0(metric, "_obs_pvalue_lower"))), NA_character_,
        ifelse(is.na(get(paste0(metric, "_obs_pvalue_upper"))), NA_character_,
               ifelse(get(paste0(metric, "_obs_pvalue_lower")) > 0.99, "a_< 0.01",
                      ifelse(get(paste0(metric, "_obs_pvalue_lower")) > 0.975, "b_< 0.025",
                             ifelse(get(paste0(metric, "_obs_pvalue_upper")) > 0.99, "e_> 0.99",
                                    ifelse(get(paste0(metric, "_obs_pvalue_upper")) > 0.975, "d_> 0.975",
                                           "c_not significant"
                                    )))))))
    return(res)
}
