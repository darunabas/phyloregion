#' Create illustrative sparse matrix
#'
#' This function visualizes a sparse matrix using vertical bands corresponding
#' to presence or absence of a species in an area.
#' @param x A matrix
#' @param col A vector of colors to represent presence or absence of a species
#' @param lwd Line width
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot.sparse
#' @importFrom graphics abline clip
#' @return Returns no value, just plot sparse matrix
#' @export
#' @method plot sparse
#' @rawNamespace export(plot.sparse)
plot.sparse <- function(x, col = c("red", "yellow"), lwd = 1, ...) {
  x <- as.matrix(x)
  #x <- x[nrow(x):1, ]
  plot(1, 1, xlab = NA, ylab = NA, ylim = c(0,1), xlim = c(0, ncol(x)),
       type = "n", xaxt = "n", yaxt = "n", ...)
  usr <- par("usr")
  x <- as.matrix(x)
  hh <- rev(seq(0, usr[4], length.out = (nrow(x) + 1)))
  states <- c(1, 0)
  suppressWarnings(invisible(lapply(seq_len(nrow(x)), function(i) {
    pos <- seq_along(x[i,])
    ind <- match(as.numeric(x[i,]), states)
    clip(usr[1], usr[2], hh[i+1], hh[i])
    abline(v = pos-0.5, col = col[ind], lwd = lwd)
    do.call("clip", as.list(usr))
  })))
}







