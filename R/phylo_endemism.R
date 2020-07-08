#' Phylogenetic Endemism
#'
#' Calculates phylogenetic endemism (sum of 'unique' branch lengths) of multiple
#' ecological samples.
#'
#' Takes a community data table and a (rooted) phylogenetic tree (with branch
#' lengths) and calculates either strict or weighted endemism in Phylogenetic
#' Diversity (PD). Strict endemism equates to the total amount of branch length
#' found only in the sample/s and is described by Faith et al. (2004) as
#' PD-endemism. Weighted endemism calculates the "spatial uniqueness" of each
#' branch in the tree by taking the reciprocal of its range, multiplying by
#' branch length and summing for all branch lengths present at a sample/site.
#' Range is calculated simply as the total number of samples/sites at which the
#' branch is present. This latter approach is described by Rosauer et al. (2009)
#' as Phylogenetic endemism.
#'
#' @param x is the community data given as a data.frame or matrix with
#' species/OTUs as columns and samples/sites as rows (like in the vegan
#' package). Columns are labeled with the names of the species/OTUs. Rows are
#' labelled with the names of the samples/sites. Data can be either abundance or
#' incidence (0/1). Column labels must match tip labels in the phylogenetic tree
#' exactly!
#' @param  phy a (rooted) phylogenetic tree (phylo) with branch lengths
#' @param weighted is a logical indicating whether weighted endemism (default)
#' or strict endemism should be calculated.
#' @return \code{phylo_endemism} returns a vector of phylogenetic endemism for
#' each sample or site.
#' @references
#' Faith, D.P., Reid, C.A.M. & Hunter, J. (2004) Integrating phylogenetic
#' diversity, complementarity, and endemism for conservation assessment.
#' \emph{Conservation Biology} \strong{18}(1): 255-261.
#'
#' Rosauer, D., Laffan, S.W., Crisp, M.D., Donnellan, C. & Cook, L.G. (2009).
#' Phylogenetic endemism: a new approach for identifying geographical
#' concentrations of evolutionary history. \emph{Molecular Ecology}
#' \strong{18}(19): 4061-4072.
#'
#' Daru, B.H., Farooq, H., Antonelli, A. & Faurby, S. (2020) Endemism
#' patterns are scale dependent. \emph{Nature Communications} \strong{11}
#' : 2115.
#' @examples
#' data(africa)
#' pe <- phylo_endemism(africa$comm, africa$phylo)
#' plot(density(pe))
#' @importFrom Matrix Diagonal crossprod
#' @importFrom phangorn Ancestors
#' @export
phylo_endemism <- function(x, phy, weighted = TRUE){
    if(length(setdiff(colnames(x), phy$tip.label)) > 0)
      stop("There are species labels in community matrix missing in the tree!")
    if(length(setdiff(phy$tip.label, colnames(x))) > 0)
        phy <- keep.tip(phy, intersect(phy$tip.label, colnames(x)))
    comm_phylo <- phylo_community(x, phy)
    weights <- comm_phylo$Matrix %*%
        Diagonal(x = 1 / colSums(comm_phylo$Matrix) )
    if (weighted == FALSE) weights[weights < 1] <- 0
    pd <- (weights %*% comm_phylo$edge.length)[,1]
    pd <- pd[row.names(x)]
    return(pd)
}
