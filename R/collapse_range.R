#' Collapse nodes and ranges based on divergence times
#'
#' This function collapses nodes and geographic ranges based on species'
#' divergence times at various time depths.
#'
#' @param x A community matrix or data frame.
#' @param tree A phylogenetic tree.
#' @param n Time depth to slice the phylogenetic tree (often in millions of
#' years for dated trees).
#' @param species If \code{format =} \dQuote{long} (the default),
#' the column with the species name.
#' @param grids The column with the sites or grids if \code{format =} \dQuote{long}.
#' @param format Format of the community composition data:
#' \dQuote{long} or \dQuote{wide} with species as columns and sites as rows.
#' @rdname collapse_range
#' @importFrom ape keep.tip
#'
#' @return Two community data frames: the \code{collapsed community data} and
#' \code{original community data}
#'
#' @references
#' Daru, B.H., Farooq, H., Antonelli, A. & Faurby, S. Endemism patterns are
#' scale dependent.
#' @examples
#' library(ape)
#' tr1 <- read.tree(text ="(((a:2,(b:1,c:1):1):1,d:3):1,e:4);")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 1,0,1,1,0,1,
#'                 0,0,0,1,1,0), 6, 5,
#'               dimnames=list(paste0("g",1:6), tr1$tip.label))
#'
#' collapse_range(com, tr1, n=1)
#' @export collapse_range
collapse_range <- function(x, tree, n, species="species", grids="grids",
                           format="wide"){
  if (format == "wide") {
    x <- data.frame(expand.grid(dimnames(provideDimnames(x)))[1:2],
                    as.vector(as.matrix(x)))
    x <- x[(x[, 3] > 0) & !is.na(x[, 3]), ]
    x <- x[sort.list(x[, 1]), ]
    x <- data.frame(grids = x[, 1], species = x[, 2])
    x <- as.data.frame(x)}
  else if (format == "long") {
    x <- as.data.frame(x)
    x <- x[, c(grids, species)]
    x <- x[, c(grids, species)]
    names(x) <- c("grids", "species")}
  ind <- intersect(x$species, tree$tip.label)
  if(length(ind)>0){
    subphy <- keep.tip(tree, ind)
    submat <- subset(x, x$species %in% ind)
    d <- get_clades(subphy, n)

    # SPECIES THAT CAN BE LUMPED
    m <- do.call("rbind", lapply(d, function(y) cbind(species=y,
                                                      rep=parent.frame()$i[],
                                                      grouping="lumper",
                                                      ntax=length(y))))
    r <- data.frame(m[m[, "ntax"] > 1,])

    # MASTER_DISTINCT -- SPECIES WHICH ARE DISTINCT AND CAN'T BE LUMPED
    w <- subset(submat, !c(submat$species %in% as.vector(r$species)))

    S <- unique(as.character(r$rep))

    out <- NULL
    for (i in seq_along(S)) {
      v <- subset(r, r$rep==S[i])
      v <- subset(submat, c(submat$species %in% as.vector(v$species)))
      v <- cbind(grids=v[1], species=v$species[1])
      out <- rbind(out, v)
    }
    rr <- unique(data.frame(out))
    z <- rbind(w, rr)
    return(list(collapsed_com=z, original_com=submat))
  } else {
    stop("Taxa names in community data do not match names in tree")
  }
}
