#' Create a subtree with largest overlap from a species list.
#'
#' phylobuilder creates a subtree with largest overlap from a species list.
#' If species in the species list are not already in the tip label, species will
#' be added at the most recent common ancestor at the genus or family level when
#' possible.
#'
#' @param species A vector or matrix containing a species list
#' @param tree a phylogenetic tree (object of class phylo)
#' @seealso \code{\link[phangorn]{add.tips}}, \code{\link[ape]{label2table}},
#' \code{\link[ape]{stripLabel}}
#' @keywords bioregion
#' @importFrom ape drop.tip stripLabel label2table getMRCA
#' @importFrom parallel mclapply
#' @importFrom phangorn add.tips
#' @examples
#' data("africa")
#' phylobuilder(species = colnames(africa$comm), tree = africa$phylo)
#' @export
phylobuilder <- function(species, tree) {  # , sep=NULL,
  taxonomy <- c("subspecies", "species", "genus", "family")
  #  ,"suborder", "order", "class", "phylum", "kingdom",  ) "subspecies"
  # branching phylomatic_names subspecies
  if (is.matrix(species) | is.data.frame(species)) species_list <- species
  if (is.factor(species)) species <- as.character(species)
  if (is.vector(species)) {
    species <- unique(species)
    species <- stripLabel(species, subsp = FALSE)
    species_list <- label2table(species)
    species_list$species <- species
  }

  if (is.matrix(species_list) | is.data.frame(species_list)) {
    colnames(species_list) <- tolower(colnames(species_list))
    nam <- intersect(taxonomy, colnames(species_list))
    species_list <- species_list[, nam]
  }

  species_list[vapply(species_list, is.factor, FALSE)] <-
    lapply(species_list[vapply(species_list, is.factor, FALSE)], as.character)

  if (any(duplicated(species_list$species))) {
    warning("Duplicated species detected and removed.")
    species_list <- species_list[!duplicated(species_list$species), ]
  }

  species_list$species <- gsub(" ", "_", species_list$species)
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  #  tree$tip.label <- stripLabel(tree$tip.label)

  species_list$species <- gsub("(^[[:alpha:]])", "\\U\\1", species_list$species,
    perl = TRUE)
  species_list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", species_list$genus,
    perl = TRUE)
  if (!is.null(species_list$family)) species_list$family <-
    gsub("(^[[:alpha:]])", "\\U\\1", species_list$family, perl = TRUE)

  T1 <- tree
  add_tip <- which(is.na(match(species_list$species, tree$tip.label)))

  species_list$add <- rep("tree", nrow(species_list))
  species_list$add[add_tip] <- "missing"

  if (length(add_tip) > 0) {

    unique_genera <- unique(species_list$genus[add_tip])
    tree_genera <- as.character(label2table(T1$tip)$"genus")

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
      tips2add <- species_list$species[add_tip[ll == 1]]
    }
    if (any(ll > 1)) {
      tips2add <- c(tips2add, species_list$species[add_tip[ll > 1]])
      fun2 <- function(x, tree) ifelse(length(x) > 1, getMRCA(tree, x), 0)
      num <- mclapply(nn_unique,  fun2, tree = T1)
      num <- (unlist(num)[tmp])[ll > 1]
      where2add <- c(where2add, unlist(num))
    }
    if (!is.null(where2add)) T1 <- add.tips(T1, tips2add, where = where2add)
    species_list$add[add_tip[ll > 0]] <- "genus"
    add_tip <- add_tip[ll == 0]

    if (length(add_tip) > 0) {
      if ("family" %in% colnames(species_list)) {
        ind <- add_tip [species_list$family[add_tip] %in% T1$node.label]
        if (length(ind) > 0) {
          num <- match(species_list$family[ind], T1$node.label)
          +length(T1$tip.label)
          T1 <- add.tips(T1, tips = species_list$species[ind], num)
          species_list$add[ind] <- "family"
        }
      }
    }

  }
  toDrop <- setdiff(T1$tip.label,
    species_list$species[species_list$add != "missing"])
  T1 <- drop.tip(T1, tip = toDrop)
  attr(T1, "species_list") <- species_list
  T1
}
