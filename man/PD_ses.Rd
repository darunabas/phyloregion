% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PD_ses.R
\name{PD_ses}
\alias{PD_ses}
\title{Phylogenetic diversity standardized for species richness}
\usage{
PD_ses(
  x,
  phy,
  model = c("tipshuffle", "rowwise", "colwise"),
  reps = 10,
  metric = "pd",
  ...
)
}
\arguments{
\item{x}{a (sparse) community matrix, i.e. an object of class matrix or
Matrix.}

\item{phy}{a phylogenetic tree (object of class phylo).}

\item{model}{The null model for separating patterns from processes and
for contrasting against alternative hypotheses. Available null models
include:
\itemize{
\item \dQuote{tipshuffle}: shuffles tip labels multiple times.
\item \dQuote{rowwise}: shuffles sites (i.e., varying richness) and
keeping species occurrence frequency constant.
\item \dQuote{colwise}: shuffles species occurrence frequency and
keeping site richness constant.}}

\item{reps}{Number of replications.}

\item{metric}{The phylodiversity measure to compute.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
A data frame of results for each community or grid cell
\itemize{
  \item grids: Site identity
  \item richness: Number of taxa in community
  \item pd_obs: Observed PD in community
  \item pd_rand.mean: Mean PD in null communities
  \item pd_rand.sd: Standard deviation of PD in null communities
  \item pd_obs.rank: Rank of observed PD vs. null communities
  \item pd_obs.z: Standardized effect size of PD vs. null communities
  \eqn{= (pd_obs - pd_rand.mean) / pd_rand_sd}
  \item pvalue: P-value (quantile) of observed PD vs. null communities
  \eqn{= mpd_obs_rank / iter + 1}
  \item reps: Number of replicates
  \item p_obs_c_lower: Number of times observed value < random value
  \item p_obs_c_upper: Number of times observed value > random value
  \item p_obs_p_lower: Percentage of times observed value < random value
  \item p_obs_p_upper: Percentage of times observed value > random value
  \item p_obs_q: Number of the non-NA random values used for comparison
}
}
\description{
This function computes the standard effect size of PD by correcting for
changes in species richness. The novelty of this function is its ability
to utilize sparse community matrix making it possible to
efficiently randomize very large community matrices spanning thousands of
taxa and sites.
}
\examples{
library(ape)
library(Matrix)
tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
  c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
  dimnames = list(paste0("g", 1:6), tree$tip.label))

PD_ses(com, tree, model="rowwise")

}
\references{
Proches, S., Wilson, J.R.U. & Cowling, R.M. (2006) How much evolutionary
history in a 10 x 10m plot? \emph{Proceedings of Royal Society B}
\strong{273}: 1143-1148.
}
