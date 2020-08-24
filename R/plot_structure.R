#plot.barnabas <- function(obj, type="pies") {
#    if(type=="structure") {
#        plot_structure()
#    }
#}

#' Visualize biogeographic patterns using structure plots
#'
#' @param omega a matrix of phyloregion of probabilities of each species
#' @param shp a polygon shapefile of grid cells with a column labeled
#' \dQuote{grids}.
#' @param col List of colors for the pies.
#' @param by Specification of the column used for grouping.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname plot_structure
#' @importFrom graphics Axis barplot layout
#' @importFrom utils modifyList
#' @return Returns no value, just visualize structure plot!
#' @export
plot_structure <- function(omega, shp, by = NULL, col=hcl.colors(12), ...) {

    #shp <- obj$shp
    #omega <- obj$omega
    index <- intersect(shp$grids, row.names(omega))

    s1 <- subset(shp, shp$grids %in% index)

    if (!is.null(by)) {
        s1 <- s1[, c("grids", by)]
        names(s1) <- c("grids", "by")
        sx <- split(s1, f=s1$by)
        vv <- as.numeric(matrix(lapply(sx, function(x) nrow(x)/nrow(s1))))

        nf <- layout(matrix(1:length(sx), ncol=1), heights=vv)

        suppressWarnings(invisible(lapply(sx, function(x) {
            ind1 <- intersect(x$grids, row.names(omega))
            y1 <- omega[ind1,]
            z1 <- t(y1[nrow(y1):1, ])
            par(mgp=c(3,0,0))
            par(mar=c(0, 5, 0, 1))
            barplot(z1, horiz=TRUE, col=col,
                    border=NA, space=0, las=1,
                    axes=FALSE, ylab=x$by[1], cex.lab=0.75)

        })))
    } else {
        ind1 <- intersect(s1$grids, row.names(omega))
        y1 <- omega[ind1,]
        z1 <- t(y1[nrow(y1):1, ])
        barplot(z1, horiz=TRUE, col=col,
                border=NA, space=0, las=1,
                axes=FALSE, cex.lab=0.75)
    }
    Axis(side=1, labels=TRUE, lwd=1, las=1, cex.axis = 1, line=-0.25)

}


