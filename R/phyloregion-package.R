#' Plants of southern Africa
#'
#' This dataset consists of a dated phylogeny of the woody plant
#' species of southern Africa along with their geographical distributions.
#' The dataset comes from a study that maps tree diversity hotspots in southern Africa
#' (Daru et al. 2015). The study mapped five types of diversity hotspots including
#' species richness (SR), phylogenetic diversity (PD), phylogenetic endemism
#' (PE), species weighted endemism (CWE), and evolutionary distinctiveness and
#' global endangerment (EDGE). The results revealed large spatial incongruence
#' between biodiversity indices, resulting in unequal representation of PD,
#' SR, PE, CWE and EDGE in hotspots and currently protected areas, suggesting
#' that an integrative approach which considers multiple facets of biodiversity
#' is needed to maximise the conservation of tree diversity in southern Africa.
#' Specifically for this package, we arranged the dataset into four components:
#' \dQuote{comm}, \dQuote{polys}, \dQuote{phylo}, \dQuote{mat}, \dQuote{IUCN}.
#'
#' @section Details:
#' \itemize{
#'   \item comm: This a sparse community composition matrix of each species
#'   presences/absences within 50 × 50 km grid cells. A sparse matrix is a
#'   matrix with a high proportion of zero entries (Duff 1977), of which only
#'   the non-zero entries are stored and used for downstream analysis.
#'   \item polys: These are the grid cells covering the geographic extent of study area.
#'   These can be created using the function \code{\link{fishnet}}. The polys object is
#'   of class \code{SpatialPolygonsDataFrame} and has a column labeled \dQuote{grids}, with
#'   the grid identities.
#'    \item phylo: This corresponds to the phylogenetic tree which was
#'    estimated using Bayesian analysis of 1,400 species and 1,633 bp
#'    of chloroplast DNA sequences derived from a combination of  \emph{matK}
#'    and  \emph{rbcLa}, assuming an uncorrelated relaxed molecular clock model,
#'    using the program BEAST v.1.7.5 (Drummond & Rambaut, 2007). Branch
#'    lengths were calibrated in millions of years using a Bayesian MCMC
#'    approach by enforcing topological constraints assuming APG III backbone
#'    from Phylomatic v.3 (Webb & Donoghue, 2005) and 18 fossil calibration
#'    points from Bell et al. (2010).
#'    \item mat: This is a distance matrix of phylogenetic beta diversity between
#'    all grid cells at the 50 × 50 km scale.
#'    \item IUCN: This is a dataframe of IUCN conservation status of each woody
#'    species (LC, NT, VU, EN, CR). This is useful for analysis of Evolutionary
#'    Distinctiveness and Global Endangerment using the function \code{\link{EDGE}}.
#' }
#' @references
#'
#' Bell, C.D., Soltis, D.E., & Soltis, P.S. (2010). The age and diversification
#' of the angiosperms re-revisited. \emph{American Journal of Botany} \strong{97},
#' 1296–1303.
#'
#' Daru, B.H., Van der Bank, M. & Davies, T.J. (2015) Spatial incongruence
#' among hotspots and complementary areas of tree diversity in southern Africa.
#' \emph{Diversity and Distributions} \strong{21}, 769-780.
#'
#' Drummond, A.J., & Rambaut, A. (2007). BEAST: Bayesian evolutionary analysis
#' by sampling trees. \emph{BMC Evolutionary Biology} \strong{7}, 214.
#'
#' Duff, I.S. (1977). A survey of sparse matrix research. \emph{Proceedings
#' of the IEEE} \strong{65}, 500–535.
#'
#' Webb, C.O., & Donoghue, M.J. (2005). Phylomatic: Tree assembly for applied
#' phylogenetics. \emph{Molecular Ecology Notes} \strong{5}, 181–183.
#'

#' @name africa
#' @docType data
#' @keywords datasets
#' @examples
#' data(africa)
#' names(africa)
#' \donttest{
#' library(raster)
#' library(ape)
#' plot(africa$polys)
#' plot(africa$phylo)
#' }
NULL


#' Biogeographic regionalization and spatial conservation
#'
#' This document describes the \code{phyloregion} package for the R software.
#' \code{phyloregion} is a computational infrastructure for biogeographic regionalization
#' (the classification of geographical areas in terms of their biotas) and
#' spatial conservation in the R scientific computing environment. Previous
#' analyses of biogeographical regionalization were either focused on smaller
#' datasets or slower particularly when the number of species or geographic scale
#' is very large. With macroecological datasets of ever increasing size and complexity,
#' \code{phyloregion} offers the possibility of handling and executing large scale
#' biogeographic regionalization efficiently and with extreme speed. It also
#' allows fast and efficient for analysis of more standard conservation measures
#' such as phylogenetic diversity, phylogenetic endemism, evolutionary distinctiveness
#' and global endangerment. \code{phyloregion} can run on any operating system
#' (Mac, Linux, Windows or even high performance computing cluster) with
#' R 3.6.0 (or higher) installed.
#'
#' @section How to cite \code{phyloregion}:
#'
#' The original implementation of phyloregion is described in:
#'
#' \itemize{
#'   \item Daru B.H., Karunarathne, P. & Schliep, K. (2020) phyloregion:
#'   R package for biogeographic regionalization and spatial
#'   conservation. \emph{bioRxiv} 2020.02.12.945691 doi: 10.1101/2020.02.12.945691
#' }
#'
#' It is based on the method described in:
#'
#' \itemize{
#'   \item Daru, B.H., Farooq, H., Antonelli, A. & Faurby, S. (2020) Endemism
#'   patterns are scale dependent. \emph{Coming soon}.
#' }
#'
#' The original conceptual is described in:
#'
#' \itemize{
#'   \item Daru, B.H., Elliott, T.L., Park, D.S. & Davies, T.J. (2017)
#'   Understanding the processes underpinning patterns of phylogenetic
#'   regionalization. \emph{Trends in Ecology and Evolution} \strong{32}: 845-860.
#' }
#'
#' @section Feedback:
#'
#' If you have any questions, suggestions or issues regarding the package,
#' please add them to \href{https://github.com/darunabas/phyloregion/issues}{GitHub issues}
#'
#' @section Installation:
#'
#' \code{phyloregion} is an open-source and free package hosted on
#' \href{https://github.com/darunabas/phyloregion}{GitHub}.
#' You will need to install the \code{devtools} package. In \code{R}, type:
#'
#' \code{if (!requireNamespace("devtools", quietly = TRUE))
#'     install.packages("devtools")}
#'
#' Then:
#'
#' \code{devtools::install_github("darunabas/phyloregion")}
#'
#' Load the phyloregion package:
#'
#' \code{library(phyloregion)}
#'
#' @section Acknowledgments:
#' Barnabas Daru thanks Texas A&M University-Corpus Christi for financial
#' and logistic support.
#'
#'
#' @name phyloregion-package
#' @aliases phyloregion-package
#' @docType package
#' @author \href{https://barnabasdaru.com/}{Barnabas H. Daru}, Piyal Karunarathne,
#' \href{https://kschliep.netlify.com/}{Klaus Schliep}
#'
#' @keywords package
NULL


