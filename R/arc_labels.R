getYmult <- function ()
{
  if (dev.cur() == 1) {
    warning("No graphics device open.")
    ymult <- 1
  }
  else {
    xyasp <- par("pin")
    xycr <- diff(par("usr"))[c(1, 3)]
    ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
  }
  return(ymult)
}


draw_arc <- function (x = 1, y = NULL, radius = 1, angle1 = deg1 * pi/180,
                      angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 0.05,
                      col = NA, lwd = NA, ...)
{
  if (all(is.na(col)))
    col <- par("col")
  if (all(is.na(lwd)))
    lwd <- par("lwd")
  xylim <- par("usr")
  ymult <- getYmult()
  devunits <- dev.size("px")
  draw.arc.0 <- function(x, y, radius, angle1, angle2, n,
                         col, lwd, ...) {
    delta.angle <- (angle2 - angle1)
    if (n != as.integer(n))
      n <- as.integer(1 + delta.angle/n)
    delta.angle <- delta.angle/n
    angleS <- angle1 + seq(0, length = n) * delta.angle
    angleE <- c(angleS[-1], angle2)
    if (n > 1) {
      half.lwd.user <- (lwd/2) * (xylim[2] - xylim[1])/devunits[1]
      adj.angle <- delta.angle * half.lwd.user/(2 * (radius +
                                                      half.lwd.user))
      angleS[2:n] <- angleS[2:n] - adj.angle
      angleE[1:(n - 1)] <- angleE[1:(n - 1)] + adj.angle
    }
    p1x <- x + radius * cos(angleS)
    p1y <- y + radius * sin(angleS) * ymult
    p2x <- x + radius * cos(angleE)
    p2y <- y + radius * sin(angleE) * ymult
    segments(p1x, p1y, p2x, p2y, col = col, lwd = lwd, ...)
  }
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  a1 <- pmin(angle1, angle2)
  a2 <- pmax(angle1, angle2)
  angle1 <- a1
  angle2 <- a2
  args <- data.frame(x, y, radius, angle1, angle2, n, col,
                     lwd, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(args))) do.call("draw.arc.0", c(args[i,
  ], ...))
  invisible(args)
}

arctext <- function (x, center = c(0, 0), radius = 1, start = NULL,
                     middle = pi/2, end = NULL, stretch = 1, clockwise = TRUE,
                     cex = NULL, ...)
{
  oldcex <- par("cex")
  if (is.null(cex))
    cex <- oldcex
  par(cex = cex)
  xvec <- strsplit(x, "")[[1]]
  lenx <- length(xvec)
  xwidths <- stretch * strwidth(xvec)
  charangles <- xwidths/radius
  changrang <- range(charangles)
  charangles[charangles < changrang[2]/2] <- changrang[2]/2
  if (!is.null(end)) {
    if (clockwise)
      start <- end + sum(charangles)
    else start <- end - sum(charangles)
  }
  if (is.null(start)) {
    if (clockwise)
      start <- middle + sum(charangles)/2
    else start <- middle - sum(charangles)/2
  }
  if (clockwise) {
    charstart <- c(start, start - cumsum(charangles)[-lenx])
    charpos <- charstart - charangles/2
  }
  else {
    charstart <- c(start, start + cumsum(charangles)[-lenx])
    charpos <- charstart + charangles/2
  }
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
  for (xchar in 1:lenx) {
    srt <- 180 * charpos[xchar]/pi - 90
    text(center[1] + radius * cos(charpos[xchar]), center[2] +
           radius * sin(charpos[xchar]) * ymult, xvec[xchar],
         adj = c(0.5, 0.5), srt = srt + 180 * (!clockwise),
         ...)
  }
  par(cex = oldcex)
}

#' Add arc labels to plotted phylogeny
#'
#' @param phy An object of class phylo.
#' @param tips A character vector (or a list) with names of the
#' tips that belong to the clade or group. If multiple groups are
#' to be plotted, tips must be given in the form of a list.
#' @param text Desired clade label.
#' @param plot_singletons Logical. If TRUE (default), adds arcs
#' (and labels) to single tip lineages. If FALSE, no arc or
#' labels will be plotted over that tip..
#' @param ln.offset Line offset (as a function of total tree height)
#' @param lab.offset Label offset.
#' @param cex Character expansion
#' @param orientation Orientation of the text. Can be \dQuote{vertical},
#'  \dQuote{horizontal}, or \dQuote{curved}.
#' @param \dots Further arguments passed to or from other methods.
#' @export
#' @importFrom ape .PlotPhyloEnv
#' @importFrom grDevices dev.size dev.cur
arc_labels <- function(phy, tips, ...) {
  UseMethod("arc_labels", tips)
}

#' @return \code{NULL}
#'
#' @rdname arc_labels
#' @method arc_labels default
#' @examples
#' \donttest{
#' old.par <- par(no.readonly = TRUE)
#' require(ape)
#' data(africa)
#' par(mai=rep(0,4))
#' plot(africa$phylo, type = "fan", show.tip.label=FALSE,
#'      open.angle = 180, edge.width=0.5)
#'
#' y <- data.frame(species=africa$phylo$tip.label)
#' y$genus <- gsub("_.*", "\\1", y$species)
#'
#' fx <- split(y, f=y$genus)
#'
#' suppressWarnings(invisible(lapply(fx, function(x) {
#'   y <- seq(from = 1.03, to = 1.09, by = ((1.09 - 1.03)/(length(fx) - 1)))
#'   z <- sample(y, 1, replace = FALSE, prob = NULL)
#'   if(nrow(x) > 10L) arc_labels(phy = africa$phylo, tips=x$species,
#'                             text=as.character(unique(x$genus)),
#'                             orientation = "curved", cex=0.5,
#'                             lab.offset = z)
#' })))
#' par(old.par)
#' }
#' @export
arc_labels.default <- function (phy = NULL, tips, text, plot_singletons = TRUE,
                               ln.offset = 1.02, lab.offset = 1.06,
                               cex = 1, orientation = "horizontal", ...){

  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (obj$type != "fan")
    stop("method works only for type=\"fan\"")
  h <- max(sqrt(obj$xx^2 + obj$yy^2))
  if (methods::hasArg(lwd))
    lwd <- list(...)$lwd
  else lwd <- graphics::par()$lwd
  if (methods::hasArg(col))
    col <- list(...)$col
  else col <- graphics::par()$col
  if (methods::hasArg(lend))
    lend <- list(...)$lend
  else lend <- graphics::par()$lend
  if (methods::hasArg(clockwise))
    clockwise <- list(...)$clockwise
  else clockwise <- TRUE
  if (methods::hasArg(n))
    n <- list(...)$n
  else n <- 0.05
  if (is.null(phy)) {
    message("phy argumengt is NULL; tip labels are being replaced by numbers.")
    phy <- list(edge = obj$edge, tip.label = 1:obj$Ntip,
                Nnode = obj$Nnode)
    class(phy) <- "phylo"
  }
  if(is.character(tips)){
    d <- which(phy$tip.label%in%tips)
  }
  d <- sort(d[d <= ape::Ntip(phy)])
  breaks <- c(0, which(diff(d) != 1), length(d))
  breaks <- lapply(seq(length(breaks) - 1),
                   function(i) d[(breaks[i] + 1):breaks[i+1]])
  for (i in seq(length(breaks))){
    if(!plot_singletons & length(breaks[[i]]) > 1) {
      next }
    d <- breaks[[i]]
    deg <- atan(obj$yy[d]/obj$xx[d]) * 180/pi
    # print(deg)
    ii <- intersect(which(obj$yy[d] >= 0), which(obj$xx[d] <
                                                   0))
    deg[ii] <- 180 + deg[ii]
    ii <- intersect(which(obj$yy[d] < 0), which(obj$xx[d] < 0))
    deg[ii] <- 180 + deg[ii]
    ii <- intersect(which(obj$yy[d] < 0), which(obj$xx[d] >=
                                                  0))
    deg[ii] <- 360 + deg[ii]
    # ln.offset <- 1.16
    # message(paste0("deg", deg, "\n"))
    draw_arc(x = 0, y = 0, radius = ln.offset * h, deg1 = min(deg),
             deg2 = max(deg), lwd = lwd, col = col, lend = lend, n = n)
    if(!is.null(text)){
      if (orientation == "curved")
        arctext(text, radius = lab.offset * h,
                middle = stats::median(range(deg * pi/180)),
                cex = cex, clockwise = clockwise)
      else if (orientation == "horizontal") {
        if(methods::hasArg(label_degree)){
          label_degree <- list(...)$label_degree
          #print(label_degree)
          if(is.null(label_degree) | is.na(label_degree)){
            label_degree <- stats::median(deg)
          }
        } else {
          label_degree <- stats::median(deg)
        }
        x0 <- lab.offset * cos(label_degree * pi/180) * h
        y0 <- lab.offset * sin(label_degree * pi/180) * h
        graphics::text(x = x0, y = y0, label = text,
                       adj = c(if (x0 >=0) 0 else 1, if (y0 >= 0) 0 else 1),
                       offset = 0, cex = cex)
      }
    }
  }
}


# #' @return \code{NULL}
# #'
# #' @rdname arclabels
# #' @method arclabels list
# #' @export
# arclabels.list <- function (phy = NULL, text, tips, plot_singletons = TRUE,
# ln.offset = 1.02, lab.offset = 1.06,
#     cex = 1, orientation = "horizontal", ...){
