% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo_beta.R
\name{phylobeta_core}
\alias{phylobeta_core}
\alias{phylobeta}
\title{Phylogenetic beta diversity}
\usage{
phylobeta_core(x, phy)

phylobeta(x, phy, index.family = "sorensen")
}
\arguments{
\item{x}{an object of class Matrix, matrix or phyloseq}

\item{phy}{a phylogenetic tree (object of class phylo)}

\item{index.family}{family of dissimilarity indices, partial match of
"sorensen" or "jaccard".}
}
\value{
\code{phylobeta_core} returns an object of class "phylo.betapart",
see \code{\link{phylo.betapart.core}} for details. This object can be called
by \code{\link{phylo.beta.pair}} or \code{\link{phylo.beta.multi}}.

\code{phylobeta} returns a list with three phylogenetic dissimilarity
matrices. See \code{\link{phylo.beta.pair}} for details.
}
\description{
\code{phylobeta_core} computes efficiently for large community matrices and
trees the necessary quantities used by the betapart package to compute
pairwise and multiple-site phylogenetic dissimilarities.
}
\examples{
library(ape)
tree <- read.tree(text = "((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
  c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
  dimnames = list(paste0("g", 1:6), tree$tip.label))
com

pbc <- phylobeta_core(com, tree)
pb <- phylobeta(com, tree)
}
\seealso{
\code{\link{read.community}}, \code{\link{phylo.betapart.core}},
\code{\link{beta_core}}
}
\author{
Klaus Schliep
}
\keyword{phyloregion}
