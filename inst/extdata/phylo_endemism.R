.SpecRich <- function(x){
  mm1 <- data.frame(table(x$grids))
  names(mm1) <- c("grids", "SR")
  mm1
}

#' Measure the degree to which phylogenetic diversity is restricted to
#' an area.
#'
#' \code{phylo.endemism} is species richness inversely weighted
#' by species ranges.
#'
#' @param x A community matrix or data frame.
#' @rdname phylo_endemism
#' @keywords bioregion
#' @importFrom raster raster rasterToPolygons xyFromCell ncell
#' @importFrom raster values
#' @importFrom sp CRS proj4string
#' @importFrom data.table as.data.table
#'
#' @return A data frame of species traits by site.
#' \deqn{ \frac{X^2}{Y} }{ X2 / Y }
#'
#' @references
#' Crisp, M.D., Laffan, S., Linder, H.P. & Monro, A. (2001) Endemism in the Australian flora.
#' \emph{Journal of Biogeography} \strong{28}: 183–198.
#'
#' Laffan, S.W., & Crisp, M.D. (2003) Assessing endemism at multiple spatial scales,
#' with an example from the Australian vascular flora. \emph{Journal of Biogeography} \strong{30}: 511–520.
#' @examples
#' require(raster)
#' data(africa)
#'
#' Endm <- weighted.endemism(africa$comm)
#' m <- merge(africa$polys, Endm, by="grids")
#' m <- m[!is.na(m@data$WE),]
#'
#' plot.swatch(m, values = m$WE, k=20)
#'
#' @export phylo.endemism

weighted.endemism <- function(x){
  tmp <- .SpecRich(x)
  index <- match(x$grids, tmp$grids)
  SR <- tmp$SR[index]
  ff <- table(x$species)
  x$WE <- SR/ff[x$species]
  tmp <- as.data.table(x)
  res <- tmp[, sum(WE), by=grids]
  res <- res[, list(grids, WE=V1)]
  res
}
phylo_endemism <- function (x, phy, weighted=T) {

# where x is a community data table (as in the vegan package) with species/OTUs as columns and samples/sites as rows. Columns are labelled with the names of the species/OTUs. Rows are labelled with the names of the samples/sites.
# where phy is a phylogenetic tree stored as a phylo object (as in the ape package) with terminal nodes labelled with names matching those of the community data table. Note that the function trims away any terminal taxa not present in the community data table, so it is not necessary to do this beforehand.
# where weighted=T is a boolean indicating whether weighted endemism or standard endemism is calculated.

library(ape)

### step 1: trimming the tree to match the community data table thus creating a "community tree" (sensu Cam Webb).

if (length(phy$tip.label) > length(x[1,])) {
	phy <- drop.tip (phy, which(!phy$tip.label %in% colnames(x))) }

# script is modified from that of Paradis (2006) "Analysis of Phylogenetics and Evolution with R", Springer.

### step 2: converting a community tree into a MRP matrix

# A MRP matrix, used in supertree methods, is where the membership of an OTU in a clade spanned by a branch is indicated by a 0 (no) or 1 (yes). Unlike supertree MRP matrices, our matrix includes terminal branches.
# the new phylo object model for ape 2.0 broke the original code. The following replaces that code.

phylomatrix <- matrix (0, length(phy$tip.label), length(phy$edge.length))
for (i in 1:length(phy$tip.label)) {
	lineage <- which (phy$edge[,2] == i)
	node <- phy$edge[lineage,1]
	while (node > length(phy$tip.label)+1) {
		branch <- which (phy$edge[,2] == node)
		lineage <- c (lineage, branch)
		node <- phy$edge[branch,1]
		}
	phylomatrix[i,lineage] = 1
	}

# this script fills a matrix of 0's with 1's corresponding to the incidence of an OTU on a particular branch.
# the code is pretty slow on large trees.

### step 3: re-ordering the OTUs to match.

phylomatrix <- phylomatrix[sort.list(phy$tip.label), ]
x <- x[ ,sort.list(colnames(x))]
x <- as.matrix (x)

# data are sorted to a common ordering standard, that is alphabetic order, so that OTUs match up.

### step 4: creating a community phylogeny matrix indicating incidence of branches per site

commphylo <- x %*% phylomatrix
commphylo <- ifelse (commphylo > 0, 1, 0)

### step 5: calculate the range (in sites) of each branch

ranges <- colSums(commphylo)

### step 6: calculate weightings per branch

weights <- t(commphylo) / ranges
if(weighted==F) weights <- ifelse(weights<1,0,1)

### step 6: calculating PD per sample

pd <- t(weights) %*% phy$edge.length

return (pd)

}



phylo.endemism (mydata, mytree, weighted = T)
