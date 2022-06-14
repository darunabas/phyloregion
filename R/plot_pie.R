legend_box <- function (x, y = NULL, maxradius, mab = 1.2, inset = 0, double = F)
{
  auto <- if (is.character(x))
    match.arg(x, c("bottomright", "bottom", "bottomleft",
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  asp <- get.asp()
  h <- mab * 2 * maxradius
  w <- h * asp
  if (double)
    h <- h * 2
  usr <- par("usr")
  inset <- rep(inset, length.out = 2)
  if (!is.na(auto)) {
    insetx <- inset[1L] * (usr[2L] - usr[1L])
    left <- switch(auto, bottomright = , topright = , right = usr[2L] -
                     w - insetx, bottomleft = , left = , topleft = usr[1L] +
                     insetx, bottom = , top = , center = (usr[1L] + usr[2L] -
                                                            w)/2)
    insety <- inset[2L] * (usr[4L] - usr[3L])
    top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] +
                    h + insety, topleft = , top = , topright = usr[4L] -
                    insety, left = , right = , center = (usr[3L] + usr[4L] +
                                                           h)/2)
  }
  else {
    left <- x - 1.2 * asp * maxradius
    top <- y + 1.2 * maxradius
  }
  return(c(left, top, left + w, top - h))
}


legend_pie <- function (x, y = NULL, z = NULL, labels, rd = NULL, bty = "n",
                        mab = 1.2, bg = NULL, inset = 0, ...)
{
  if (is.null(z))
    z <- rep(1, length.out = length(labels))
  box <- legend_box(x, y, rd, mab, inset)
  if (bty == "o")
    rect(box[1], box[2], box[3], box[4], col = bg)
  x <- (box[1] + box[3])/2
  y <- box[4] + mab * rd
  add_pie(z, x, y, labels, rd, ...)
}


hue <- function(x, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
                     random=FALSE) {
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0,
            hmax <= 360, cmax <= 180, lmax <= 100,
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            x > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1),
                                   seq(-100, 100, 5),
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin &
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    lab <- as(hcl, 'LAB')
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), x, iter.max=50)
  hex(LAB(clus$centers))
}

get.asp <- function () {
    pin <- par("pin")
    usr <- par("usr")
    asp <- (pin[2]/(usr[4] - usr[3]))/(pin[1]/(usr[2] - usr[1]))
    return(asp)
}

add_pie <- function (z, x = 0, y = 0, labels = names(z), radius = 1,
                     edges = 200, clockwise = TRUE, init.angle = 90,
                     density = NULL, angle = 45, col = NULL,
                     border = NULL, lty = NULL, label.dist = 1.1, ...) {

    if (!is.numeric(z) || any(is.na(z) | z < 0))
        stop("'z' values must be positive.")
    if (is.null(labels))
        labels <- as.character(seq_along(z))
    else labels <- as.graphicsAnnot(labels)
    z <- c(0, cumsum(z)/sum(z))
    dz <- diff(z)
    nz <- length(dz)
    asp <- get.asp()
    if (is.null(col))
        col <- if (is.null(density))
            c("#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B",
              "#9E67AB", "#CE7058", "#D77FB4")
    else par("fg")
    if (!is.null(col))
        col <- rep_len(col, nz)
    if (!is.null(border))
        border <- rep_len(border, nz)
    if (!is.null(lty))
        lty <- rep_len(lty, nz)
    angle <- rep(angle, nz)
    if (!is.null(density))
        density <- rep_len(density, nz)
    twopi <- if (clockwise)
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = asp * radius * cos(t2p) + x,
             y = radius * sin(t2p) + y)
    }
    for (i in 1L:nz) {
        n <- max(2, floor(edges * dz[i]))
        P <- t2xy(seq.int(z[i], z[i + 1], length.out = n))
        polygon(c(P$x, 0 + x), c(P$y, 0 + y), density = density[i],
                angle = angle[i], border = border[i], col = col[i],
                lty = lty[i])
        P <- t2xy(mean(z[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            text(label.dist * (P$x - x) + x,
                 label.dist * (P$y - y) + y,
                 labels[i], xpd = TRUE,
                 adj = ifelse(P$x - x < 0, 1, 0), ...)
        }
    }
}


#' Visualize biogeographic patterns using pie charts
#'
#' @param omega a matrix of phyloregion of probabilities of each species
#' @param shp a polygon shapefile of grid cells with a column labeled
#' \dQuote{grids}.
#' @param pie_control The list of control parameters to be passed into
#' the add.pie function.
#' @param legend_pie Legend for the pie plots.
#' @param r Radius of the pie legend to be displayed
#' @param legend Logical, whether to plot a legend or not.
#' @param col List of colors for the pies.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot_pie
#' @importFrom graphics polygon par legend
#' @importFrom sp coordinates
#' @importFrom methods slot
#' @importFrom colorspace LAB hex coords
#' @importFrom utils modifyList
#' @importFrom stats kmeans
#' @return Returns no value, just map color pies in geographic space!
#' @examples
#' data(africa)
#' K <- ncol(africa$omega)
#'
#' COLRS <- phyloregion:::hue(K)
#' plot_pie(africa$omega, shp = africa$polys, col=COLRS)
#' @export
plot_pie <- function (omega, shp, r = 1, col=hcl.colors(5),
                      pie_control = list(), legend = FALSE,
                      legend_pie = FALSE, ...) {

    index <- intersect(shp$grids, rownames(omega))
    s <- subset(shp, shp$grids %in% index)
    omega <- omega[index,]

    pie_control_default <- list(edges = 200, clockwise = TRUE,
                                init.angle = 90, density = NULL,
                                angle = 45, border = NA,
                                lty = NULL, label.dist = 1.1)
    pie_control <- modifyList(pie_control_default, pie_control)


    #COLRS <- hue(K, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
    #                 random=FALSE)

    plot(s, border = NA, ...)
        suppressWarnings(invisible(lapply(1:dim(omega)[1], function(r) {
        do.call(add_pie, append(list(
            z = as.integer(100 * omega[r, ]),
            x = coordinates(s[r,])[, 1],
            y = coordinates(s[r,])[, 2],
            labels = c("", "", ""),
            radius = sqrt(sapply(slot(s[r,], "polygons"),
                                 function(i) slot(i, "area")))*0.55,
            col = col), pie_control))
    })))
        if (isTRUE(legend)) {
          legend("bottomright", legend=colnames(omega), y.intersp = 0.8, bty = "n",
                 col = col, ncol = 2, pch = 19, pt.cex = 1.5, ...)
        }
        if (isTRUE(legend_pie)) {
            rl <- (sqrt(sapply(slot(s[1,], "polygons"),
                              function(i) slot(i, "area")))*r)*2
            legend_pie("left", labels=colnames(omega), rd=rl, bty="n",
                       col=col, cex=0.5, label.dist=1.3, border = NA)
        }

}


