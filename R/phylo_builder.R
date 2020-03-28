#' Create a subtree with largest overlap from a species list.
#'
#' phylobuilder creates a subtree with largest overlap from a species list.
#' If species in the species list are not already in the tip label, species will
#' be added at the most recent common ancestor at the genus or family level when
#' possible.
#'
#' @param species A vector or matrix containing a species list
#' @param tree a phylogenetic tree (object of class phylo)
#' @param extract extract the species in the list after trying to  add missing
#' labels to the tree. If FALSE phylobuilder adds only the taxa in the list.
#' @return \code{phylobuilder} returns a phylogenetic tree, i.e. an object of
#' class \code{phylo}.
#' @seealso \code{\link[phangorn]{add.tips}}, \code{\link[ape]{label2table}},
#' \code{\link[ape]{stripLabel}}
#' @keywords bioregion
#' @importFrom ape drop.tip stripLabel label2table getMRCA
#' @importFrom parallel mclapply
#' @importFrom phangorn add.tips
#' @examples
#' library(ape)
#' txt <- "(((((Panthera_leo,Panthera_pardus), Panthera_onca),(Panthera_uncia,
#'   (Panthera_tigris_altaica, Panthera_tigris_amoyensis)))Panthera)Felidae,
#'   (((((((Canis_lupus,Canis_lupus_familiaris),Canis_latrans),Canis_anthus),
#'   Canis_aureus),Lycaon_pictus),(Canis_adustus,Canis_mesomelas))Canis)
#'   Canidae)Carnivora;"
#' txt <- gsub("[[:space:]]", "", txt)
#' cats_and_dogs <- read.tree(text=txt)
#' plot(cats_and_dogs, node.depth=2, direction="downwards")
#' nodelabels(cats_and_dogs$node.label, frame="none", adj = c(0.5, 0))
#'
#' tree <- drop.tip(cats_and_dogs, c("Panthera_uncia", "Lycaon_pictus"),
#'   collapse.singles=FALSE)
#'
#' dogs <- c("Canis_lupus", "Canis_lupus_familiaris", "Canis_latrans",
#'   "Canis_anthus", "Canis_aureus", "Lycaon_pictus", "Canis_adustus",
#'   "Canis_mesomelas")
#'
#' # try to extract tree with all 'dogs'
#' t1 <- phylobuilder(dogs, tree)
#' plot(t1, direction="downwards")
#' attr(t1, "species_list")
#'
#' # providing extra information ("Family", "Order", ...) can help
#' sp <- data.frame(Order = c("Carnivora", "Carnivora", "Carnivora"),
#'   Family = c("Felidae", "Canidae", "Canidae"),
#'   Genus = c("Panthera", "Lycaon", "Vulpes"),
#'   Species = c("uncia", "pictus", "vulpes"),
#'   Common_name = c("Snow leopard", "Africa wild dog", "Red fox"))
#' sp
#' # Now we just add some species
#' t2 <- phylobuilder(sp, tree, extract=FALSE)
#' plot(t2, direction="downwards")
#' attr(t2, "species_list")
#'
#' @export
phylobuilder <- function(species, tree, extract=TRUE) {
  taxonomy <- c("subspecies", "species", "genus", "family", "suborder",
                "order", "class", "phylum", "kingdom")
  nodetaxonomy <- c("family", "suborder", "order", "class", "phylum")
  if (is.matrix(species) | is.data.frame(species)){
    species_list <- species
    species <- NULL
  }
  if (is.factor(species)) species <- as.character(species)
  if (is.vector(species)) {
    species <- unique(species)
    species <- stripLabel(species, subsp = FALSE)
    species_list <- label2table(species)
  }
  species_list[vapply(species_list, is.factor, FALSE)] <-
    lapply(species_list[vapply(species_list, is.factor, FALSE)], as.character)
  colnames(species_list) <- tolower(colnames(species_list))
  nam <- intersect(taxonomy, colnames(species_list))
  species_list <- species_list[, nam]
  if(is.null(species)){
    species <- paste(species_list$genus, species_list$species, sep="_")
    # TODO add subspecies
  }

  if (any(duplicated(species))) {
    warning("Duplicated species detected and removed.")
    species_list <- species_list[!duplicated(species), ]
    species <- species[!duplicated(species)]
  }

  T1 <- tree

  add_tip <- which(is.na(match(species, tree$tip.label)))
  add <- rep("tree", nrow(species_list))
  add[add_tip] <- "missing"

  if (length(add_tip) > 0) {

    unique_genera <- unique(species_list$genus[add_tip])
    tree_genera <- as.character(label2table(T1$tip.label)$"genus")
    fun <- function(g, g_tree) {
      tmp <- match(g_tree, g)
      l <- length(g)
      res <- vector("list", l)
      for (i in seq_along(tmp)) {
        a <- tmp[i]
        if (!is.na(a)) res[[a]] <- c(res[[a]], i)
      }
      res
    }

    nn_unique <- fun(unique_genera, tree_genera)
    ll_unique <- lengths(nn_unique)

    tmp <- match(species_list$genus[add_tip], unique_genera)
    nn <- nn_unique[tmp]
    ll <- ll_unique[tmp]

    tips2add <- NULL
    where2add <- NULL
    if (any(ll == 1)) {
      tmp_num <- unlist(nn[ll == 1])
      num <- T1$edge[match(tmp_num, T1$edge[, 2]), 1]
      where2add <- c(where2add, num)
      tips2add <- species[add_tip[ll == 1]]
    }
    if (any(ll > 1)) {
      tips2add <- c(tips2add, species[add_tip[ll > 1]])
      fun2 <- function(x, tree) ifelse(length(x) > 1, getMRCA(tree, x), 0)
      num <- mclapply(nn_unique,  fun2, tree = T1)
      num <- (unlist(num)[tmp])[ll > 1]
      where2add <- c(where2add, unlist(num))
    }
    if (!is.null(where2add)) T1 <- add.tips(T1, tips2add, where = where2add)
    add[add_tip[ll > 0]] <- "genus"
    add_tip <- add_tip[ll == 0]

    if ( (length(add_tip) > 0) && !is.null(tree$node.label) ) {
      nodetaxonomy <- c("genus", "family", "suborder", "order", "class",
                        "phylum")
      nodetaxonomy <- intersect(nodetaxonomy, colnames(species_list))
      for(i in nodetaxonomy){
        ind <- add_tip[species_list[add_tip, i] %in% T1$node.label]
        if (length(ind) > 0) {
          num <- match(species_list[ind, i], T1$node.label) + Ntip(T1)
          T1 <- add.tips(T1, tips = species[ind], num)
          add[ind] <- i
          add_tip <- setdiff(add_tip, ind)
        }
      }
    }
  }
  if(extract) T1 <- keep.tip(T1, species[add != "missing"])
  attr(T1, "species_list") <- cbind(species, added=add)
  T1
}

