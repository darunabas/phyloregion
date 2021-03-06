% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylobeta_ses.R
\name{phylobeta_ses}
\alias{phylobeta_ses}
\title{Phylogenetic beta diversity standardized for species beta diversity}
\usage{
phylobeta_ses(
  x,
  phy,
  index.family = "simpson",
  model = c("tipshuffle", "rowwise", "colwise"),
  reps = 1000,
  ...
)
}
\arguments{
\item{x}{a (sparse) community matrix, i.e., an object of class
matrix or Matrix.}

\item{phy}{a phylogenetic tree (object of class phylo).}

\item{index.family}{the family of dissimilarity indices including
\dQuote{simpson}, \dQuote{sorensen} and \dQuote{jaccard}.}

\item{model}{The null model for separating patterns from processes and
for contrasting against alternative hypotheses. Available null models
include:
\itemize{
\item \dQuote{tipshuffle}: shuffles phylogenetic tip labels multiple times.
\item \dQuote{rowwise}: shuffles sites (i.e., varying richness) and
keeping species occurrence frequency constant.
\item \dQuote{colwise}: shuffles species occurrence frequency and
keeping site richness constant.}}

\item{reps}{Number of replications.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
A data frame of results for each community or grid cell
\itemize{
  \item phylobeta_obs: Observed phylobeta in community
  \item phylobeta_rand_mean: Mean phylobeta in null communities
  \item phylobeta_rand_sd: Standard deviation of phylobeta in null communities
  \item phylobeta_obs_z: Standardized effect size of phylobeta vs.
  null communities \eqn{= (phylobeta_obs - phylobeta_rand_mean)
  / phylobeta_rand_sd}
  \item reps: Number of replicates
}
}
\description{
This function computes the standard effect size of phylogenetic beta
diversity by correcting for changes in species beta diversity. The
novelty of this function is its ability to utilize sparse community
matrix making it possible to efficiently randomize very large
community matrices spanning thousands of taxa and sites.
}
\examples{
library(ape)
library(Matrix)
tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
  c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
  dimnames = list(paste0("g", 1:6), tree$tip.label))

phylobeta_ses(com, tree, model="rowwise")

}
\references{
Proches, S., Wilson, J.R.U. & Cowling, R.M. (2006) How much evolutionary
history in a 10 x 10m plot? \emph{Proceedings of Royal Society B}
\strong{273}: 1143-1148.
}
