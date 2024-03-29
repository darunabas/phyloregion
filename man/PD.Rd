% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd.R
\name{PD}
\alias{PD}
\alias{RPD}
\title{Phylogenetic diversity}
\usage{
PD(x, phy)

RPD(x, phy)
}
\arguments{
\item{x}{a community matrix, i.e. an object of class matrix or Matrix or an
object of class phyloseq.}

\item{phy}{a phylogenetic tree (object of class phylo).}
}
\value{
a vector with the PD for all samples.
}
\description{
\code{PD} calculates Faith's (1992) phylogenetic diversity and relative
phylogenetic diversity.
}
\examples{
library(ape)
library(Matrix)
tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
com <- sparseMatrix(c(1,3,4,1,4,5,1,2,3,4,5,6,3,4,6),
  c(1,1,1,2,2,2,3,3,3,3,3,3,4,4,4),x=1,
  dimnames = list(paste0("g", 1:6), tree$tip.label))

PD(com, tree)
# Relative phylogenetic diversity
RPD(com, tree)
}
\references{
Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity.
\emph{Biological Conservation} \strong{61}: 1–10.
}
\seealso{
read.community read.tree phylobeta_core
}
\keyword{cluster}
