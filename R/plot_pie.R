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
#' @param x a matrix of phyloregion of probabilities of each species
#' @param shp.grids if specified, the polygon shapefile of grid cells
#' with a column labeled \dQuote{grids}.
#' @param pie_control The list of control parameters to be passed into
#' the add.pie function.
#' @param K Number of distinctive colors for the pies corresponding to the
#' number of clusters or regions.
#' @param legend Logical, whether to plot a legend or not.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot_structure
#' @importFrom graphics polygon par legend
#' @importFrom sp coordinates
#' @importFrom rgeos gArea
#' @importFrom colorspace LAB hex coords
#' @importFrom utils modifyList
#' @importFrom stats kmeans
#' @return Returns no value, just map color pies in geographic space!
#' @examples
#' data(africa)
#' plot_structure(africa$omega, shp.grids = africa$polys, K = africa$K)
#' @export
plot_structure <- function (x = NULL, shp.grids, K = 5,
                      pie_control = list(), legend = FALSE, ...) {

    index <- intersect(row.names(x), shp.grids$grids)
    s <- subset(shp.grids, shp.grids$grids %in% index)

    pie_control_default <- list(edges = 200, clockwise = TRUE,
                                init.angle = 90, density = NULL,
                                angle = 45, border = NA,
                                lty = NULL, label.dist = 1.1)
    pie_control <- modifyList(pie_control_default, pie_control)

    COLRS <- hue(K, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
                     random=FALSE)

    #old.par <- par(no.readonly = TRUE)
    #par(lwd = 0.001)
    plot(s, border = NA, ...)
        suppressWarnings(invisible(lapply(1:dim(x)[1], function(r) {
        do.call(add_pie, append(list(
            z = as.integer(100 * x[r, ]),
            x = coordinates(s[r,])[, 1],
            y = coordinates(s[r,])[, 2],
            labels = c("", "", ""),
            radius = sqrt(gArea(s[r,]))*0.5,
            col = COLRS), pie_control))
    })))
    #par(old.par)
        if (isTRUE(legend)) {
          legend("bottomleft", legend=colnames(x), y.intersp = 0.8,
                 fill = COLRS, bty = "n", ...)
        }

}


