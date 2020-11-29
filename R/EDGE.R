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
#' @importFrom stats complete.cases
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
#' y <- EDGE(x = africa$IUCN, phy = africa$phylo, Redlist = "IUCN", species="Species")
#' @export EDGE

EDGE <- function(x, phy, Redlist="Redlist", species="species", ...){

  ref_table <- data.frame(Redlist = c("LC", "NT", "VU", "EN",
                                      "CR" , "EW", "EX", "DD"),
                          GE = c(0.001, 0.01, 0.1, 0.67,
                                 0.999, 1, 1, NA))

  x <- x[, c(species, Redlist)] %>%
    setNames(c("species", "Redlist")) %>%
    left_join(ref_table, by = 'Redlist')

  # CALCULATING ED
  EDGE <- evol_distinct(phy, type = "fair.proportion") %>%
    as.data.frame() %>%
    rownames_to_column(var = 'species') %>%
    inner_join(x, by = 'species') %>%
    rename(ED = '.') %>%
    mutate(EDGE = log(1+ ED) + ( GE * log(2))) %>%
    pull(EDGE, name = species)


  return(EDGE)

}
