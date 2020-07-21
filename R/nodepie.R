utils::globalVariables(".PlotPhyloEnv")
#' Label phylogenetic nodes using pie
#'
#' @param pie Estimates from ancestral character reconstruction
#' @param pie_control The list of control parameters to be passed into
#' the add.pie function.
#' @param radius Radius of the pie
#' @param legend Logical, whether to add a legend or not.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname nodepie
#' @importFrom graphics par legend
#' @importFrom utils modifyList
#' @return Returns no value, just add color pies on phylogenetic nodes!
#' @examples
#' \donttest{
#' old.par <- par(no.readonly = TRUE)
#' require(ape)
#' data(africa)
#' m <- ace_multi(africa$theta, africa$phylo)
#'
#' par(mai=rep(0,4))
#' par(mfrow=c(2,1))
#' plot(m$tree, type = "fan", show.tip.label=FALSE,
#'      open.angle = 180, edge.width=0.5)
#'
#' nodepie(m$ancestral_states, radius = 2)
#'
#' plot_structure(africa$omega, shp = africa$polys, legend = TRUE)
#' par(old.par)
#' }
#'
#' @export
nodepie <- function (pie, radius = 2, pie_control = list(),
                     legend = FALSE, ...) {

    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    YY <- lastPP$yy[node]

    s <- cbind(XX, YY)
    s <- as.data.frame(s)

    pie_control_default <- list(edges = 200, clockwise = TRUE,
                                init.angle = 90, density = NULL,
                                angle = 45, border = NA,
                                lty = NULL, label.dist = 1.1)
    pie_control <- modifyList(pie_control_default, pie_control)

    K <- ncol(pie)

    COLRS <- hue(K, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
                 random=FALSE)

    suppressWarnings(invisible(lapply(1:dim(pie)[1], function(r) {
        do.call(add_pie, append(list(
            z = 100 * pie[r, ],
            x = s[r,][[1]],
            y = s[r,][[2]],
            labels = c("", "", ""),
            radius = radius,
            col = COLRS), pie_control))
    })))
    if (isTRUE(legend)) {
        legend("bottomright", legend=colnames(pie), y.intersp = 0.8, bty = "n",
               col = COLRS, ncol = 2, pch = 19, pt.cex = 1.5, ...)
    }
}

