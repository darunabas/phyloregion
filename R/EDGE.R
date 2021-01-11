#' Evolutionary Distinctiveness and Global Endangerment
#'
#' This function calculates EDGE by combining evolutionary distinctiveness
#' (ED; i.e., phylogenetic isolation of a species) with global endangerment
#' (GE) status as defined by the International Union for Conservation
#' of Nature (IUCN).
#'
#' EDGE is calculated as: \deqn{log(1+ED) + GE*log(2)}
#' where \emph{ED} represents the evolutionary distinctiveness score of each
#' species (function \code{evol_distinct}), i.e. the degree of phylogenetic
#' isolation, and combining it with \emph{GE}, global endangerment from IUCN
#' conservation threat categories. \emph{GE} is calculated as the expected
#' probability of extinction over 100 years of each taxon in the phylogeny
#' (Redding & Mooers, 2006), scaled as follows: least concern = 0.001, near
#' threatened and conservation dependent = 0.01, vulnerable = 0.1,
#' endangered = 0.67, and critically endangered = 0.999.
#'
#' @param x a data.frame
#' @param Redlist column in the data frame with the IUCN ranks: \code{LC},
#' \code{NT}, \code{VU}, \code{EN}, \code{CR}, and \code{EX}.
#' @param species data frame column specifying the taxon
#' @param phy a phylogenetic tree (object of class phylo).
#' @param \dots Further arguments passed to or from other methods.
#' @rdname EDGE
#' @return Returns a dataframe of EDGE scores
#' @importFrom stats complete.cases na.omit
#'
#' @author Barnabas H. Daru
#' @references
#' Redding, D.W., & Mooers, A.Ø. (2006) Incorporating evolutionary measures
#' into conservation prioritization. \emph{Conservation Biology}
#' \strong{20}: 1670–1678.
#'
#' Isaac, N.J., Turvey, S.T., Collen, B., Waterman, C. & Baillie, J.E.
#' (2007) Mammals on the EDGE: conservation priorities based on threat
#' and phylogeny. \emph{PLoS ONE} \strong{2}: e296.
#' @examples
#' data(africa)
#' y <- EDGE(x=africa$IUCN, phy=africa$phylo, Redlist="IUCN", species="Species")
#' @export EDGE

EDGE <- function(x, phy, Redlist="Redlist", species="species", ...){
  x <- as.data.frame(x)
  x <- x[, c(species, Redlist)]
  names(x) <- c("species", "Redlist")
  # Calculating GE
  lookup <- c(0.001, 0.01, 0.1, 0.67, 0.999, 1, 1, NA_real_)
  names(lookup) <- c("LC", "NT", "VU", "EN", "CR" , "EW", "EX", "DD")
  Redlist <- lookup[x$Redlist]
  names(Redlist) <- x$species
  Redlist <- na.omit(Redlist)
  # CALCULATING ED
  ED <- evol_distinct(phy, type = "fair.proportion", ...)
  index <- intersect(names(ED), names(Redlist))
  my_EDGE <- (log(1+ED[index]) + (Redlist[index] * log(2)))
  names(my_EDGE) <- index
  my_EDGE

}
