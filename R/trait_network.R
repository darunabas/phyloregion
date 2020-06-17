#' Map network of functional traits
#'
#' This function presents a network approach to functional beta diversity
#'
#' @param x A data frame, matrix or edgelist
#' @param d A data frame with the functional traits
#' @param k The number to bin discrete traits.
#' @rdname trait_network
#' @keywords phyloregion
#' @importFrom igraph graph_from_data_frame
#'
#' @return
#' {A graph object}
#'
#'
#' @examples
#' com <- data.frame(grids = c("v1", "v1", "v2", "v2", "v3", "v3", "v3"),
#'                  species = c("s1", "s3", "s2", "s3", "s1", "s2", "s3"))
#'
#' trait <- data.frame(species = c("s1", "s2", "s3"),
#'                    habit = c("tree", "herb", "liana"),
#'                    IUCN = c("LC", "VU", "CR"),
#'                    height = c(7.9, 0.5, 3))
#'
#' g <- trait_network(x=com, d=trait)
#' @export
trait_network <- function(x, d, k = 5) {
  x <- cbind(x, d[, !names(d) %in% "species"][match(x$species, d$species),])
  y <- x[ , !names(x) %in% c("species","grids", "IUCN")]
  out <- NULL
  for (i in 1:ncol(y)) {
    if(is(y[,i], "character")) {
      z <- y[,i]
      z <- cbind(1/(table(z))[z])
      colnames(z) <- names(y[i])
      out <- cbind(out, z)
    } else {
      a <- rank(y[,i], ties.method = "first")
      z <- cut(a, quantile(a, probs = seq(0, 1, 1/k), na.rm = TRUE),
                         include.lowest = TRUE, labels = FALSE)
      z <- cbind(1/(table(z))[z])
      colnames(z) <- names(y[i])
      out <- cbind(out, z)
    }
  }
  g <- cbind(x[ , names(x) %in% c("species","grids")], out)
  igraph::graph_from_data_frame(g, directed = FALSE)
}



