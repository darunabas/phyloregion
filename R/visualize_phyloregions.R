color_key <- function(y, cols, vals, lab="ED", leg=5, lwd=15, pos="bottomright") {
    if(pos=="bottomright"){
        a=y@bbox[3]+0.5
        b=y@bbox[2]
    }
    if(pos=="topleft"){
        a=y@bbox[1]-0.5
        b=y@bbox[4]-leg
    }
    if(pos=="bottomleft"){
        a=y@bbox[1]-0.5
        b=y@bbox[2]
    }
    if(pos=="topright"){
        a=y@bbox[3]+0.5
        b=y@bbox[4]-leg
    }
    X <- cbind(rep(a, length(cols)), rep(a, length(cols)))
    Y <- b + cbind(0:(length(cols) - 1)/length(cols),
                   1:length(cols)/length(cols)) * (leg)
    for (i in 1:length(cols)) lines(X[i, ], Y[i, ], col = cols[i],
                                    lwd = lwd, lend = 2)
    text(x = a, y = b, round(min(vals), 3), pos = 4, cex = 0.7) # lim texts
    text(x = a, y = b+(leg/2), round(median(vals), 3), pos = 4, cex = 0.7)
    text(x = a, y = b+leg, round(max(vals), 3), pos = 4, cex = 0.7)
    text(x = a, y = b+leg, lab, pos = 3, cex = 1)
}

#' Visualize biogeographic patterns
#'
#' @param x an object of class phyloregion from \code{ed_phyloregion}
#' @param palette name of the palette to generate colors from.
#' The name is matched to the list of available color palettes from
#' the \code{hcl.colors} function in the \code{grDevices} package.
#' @param key_label label for the color key
#' @param legend logical indicating whether to add a legend to the map.
#' @param pos location to position the legend such as \dQuote{bottomright},
#' \dQuote{bottomleft}, \dQuote{topleft}, and \dQuote{topright}.
#' @param leg parameter of the color key.
#' @param lwd parameter of the color key.
#' @param \dots arguments passed among methods.
#' @rdname plot_evoldistinct
#' @export
plot_evoldistinct <- function(x, palette = "YlOrBr", pos = "bottomleft",
                              key_label = "", legend = TRUE, leg =5,
                              lwd = 15, ...){
    if (!inherits(x, "phyloregion"))
        stop("object \"x\" is not of class \"phyloregion\"")
    m1 <- x$evol_distinct
    k <- nrow(m1)
    COLRS <- hcl.colors(k, palette, rev = TRUE, ...)
    y <- choropleth(m1, m1$ED, k) #, style = style
    plot(y, col = COLRS[y$values], ...)
    text(y, labels = as.character(y@data$cluster), ...)
    if(legend){
        color_key(y, COLRS, vals = m1$ED, leg = leg, lwd = lwd, pos = pos, lab = key_label)
    }
}


#' @rdname plot_evoldistinct
#' @importFrom raster text
#' @export
plot_phyloregion <- function(x, ...){
    if (!inherits(x, "phyloregion"))
        stop("object \"x\" is not of class \"phyloregion\"")
    y <- x$evol_distinct
    plot(y, col=y$COLOURS, ...)
    raster::text(y, labels=as.character(y@data$cluster), ...)
}

#' @rdname plot_evoldistinct
#' @export
plot_NMDS <- function(x, ...){
    if (!inherits(x, "phyloregion"))
        stop("object \"x\" is not of class \"phyloregion\"")
    c1 <- x$NMDS
    plot(c1$points, bg=hexcols(c1), pch=21, ...)
}

#' @rdname plot_evoldistinct
#' @examples
#' data(africa)
#' tree <- africa$phylo
#' x <- sampl2sparse(africa$comm)
#'
#' subphy <- match_phylo_comm(tree, x)$phy
#' submat <- match_phylo_comm(tree, x)$com
#'
#' pbc <- phylobeta(submat, subphy)
#' y <- ed_phyloregion(pbc, shp=africa$polys)
#'
#' plot_NMDS(y, cex=6)
#' text_NMDS(y, cex=2)
#' plot_phyloregion(y, cex=1)
#' plot_evoldistinct(y)
#' @export
text_NMDS <- function(x, ...){
    if (!inherits(x, "phyloregion"))
        stop("object \"x\" is not of class \"phyloregion\"")
    c1 <- x$NMDS
    text(c1$points, as.character(1:c1$nobj), ...)
}

