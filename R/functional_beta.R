.mtchtext <- function(x) {
    x <- as.data.frame(x)
    nat <- colnames(x)
    SP <- paste(c("\\bspecies\\b", "\\bbinomial\\b", "\\bbinomil\\b",
                  "\\btaxon\\b"), collapse = "|")
    species <- nat[grepl(SP, nat, ignore.case = TRUE)]
    others <- setdiff(nat, species)
    x <- x[, c(species, others)]
    names(x)[1] <- "species"
    message(paste0("Assuming `", paste(others, collapse = "`, `"),
                   "` as functional traits."))
    return(x)
}


lbow <- function(varpc, low = .08, max.pc = .9) {
    ee <- varpc/sum(varpc) # ensure sums to 1
    #print(round(log(ee),3))
    while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
    lowie <- (ee<low) ; highie <- ee>low/8
    low.ones <- which(lowie & highie)
    others <- length(which(!lowie))
    if(length(low.ones)>0) {
        if(length(low.ones)==1) {
            elbow <- low.ones
        } else {
            set <- ee[low.ones]
            pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
            infz <- is.infinite(pc.drops)
            #print(pc.drops)
            elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
        }
    } else {
        # if there are no small objective function, just choose the elbow as the second last
        cat("no objective function `E` were significantly smaller than the previous\n")
        elbow <- length(ee)
    }
    if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
        elbow <- which(cumsum(ee)>max.pc)[1]-1
    }
    if(elbow<1) {
        warning("elbow calculation failed, return zero")
        return(0)
    }
    #res <- varpc[[elbow]]
    #names(res) <- NULL
    return(list(k=as.numeric(names(varpc[elbow])), value=varpc[[elbow]]))
}


#' Functional beta diversity for mixed-type functional traits
#'
#' Computes turnover of functional diversity using k-prototypes clustering
#' algorithm tailored for mixed-type functional traits (numeric and
#' categorical) to generate an integer vector of cluster assignments. The
#' ranges of each species in a cluster are collapsed to generate a new
#' community matrix based on the presence or absence of cluster membership
#' in a grid cell. A grade of membership model or beta diversity is then fitted
#' to the new reduced community matrix for further analysis.
#'
#' @param x A dataframe or sparse community matrix of species occurrences.
#' @param trait A data frame with the first column labeled \dQuote{species}
#' containing the taxonomic groups to be evaluated whereas the remaining
#' columns contain the various functional traits. The variables should be
#' mixed-type combining numeric and categorical variables.
#' @param bin The desired number of clusters or bins. If \code{elbow=TRUE},
#' the optimal number of clusters is determined by running the analysis
#' multiple times varying from 2 to bin.
#' @param na.rm Logical, whether NA values should be removed prior to
#' computation
#' @param quick_elbow Quickly estimate the 'elbow' of a scree plot to determine
#' the optimal number of clusters.
#' @param abundance Logical, whether the reduced matrix should be returned as
#' presence or absence of cluster representation or as abundances of cluster
#' memberships
#' @param \dots Further arguments passed to or from other methods.
#' @rdname functional_beta
#' @importFrom clustMixType kproto
#' @importFrom Matrix Matrix
#' @importFrom utils tail
#' @return A list with three dissimilarity matrices capturing: (i) turnover
#' (replacement), (ii) nestedness-resultant component, and (iii) total
#' dissimilarity (i.e. the sum of both components).
#'
#' For index.family="sorensen"
#' the three matrices are:
#' \itemize{
#'   \item \code{beta.sim} A distance object, dissimilarity matrix accounting
#'   for spatial turnover (replacement), measured as Simpson pair-wise
#'   dissimilarity.
#'   \item \code{beta.sne} \code{dist} object, dissimilarity matrix accounting for
#'   nestedness-resultant dissimilarity, measured as the nestedness-fraction
#'   of Sorensen pair-wise dissimilarity
#'   \item \code{beta.sor} \code{dist} object, dissimilarity matrix accounting
#'   for total dissimilarity, measured as Sorensen pair-wise dissimilarity
#'   (a monotonic transformation of beta diversity)
#' }
#' For index.family="jaccard" the three matrices are:
#' \itemize{
#'   \item \code{beta.jtu} A distance object, dissimilarity matrix accounting
#'   for spatial turnover, measured as the turnover-fraction of Jaccard
#'   pair-wise dissimilarity
#'   \item \code{beta.jne} \code{dist} object, dissimilarity matrix accounting
#'   for nestedness-resultant dissimilarity, measured as the
#'   nestedness-fraction of Jaccard pair-wise dissimilarity
#'   \item \code{beta.jac} \code{dist} object, dissimilarity matrix accounting
#'   for beta diversity, measured as Jaccard pair-wise dissimilarity (a
#'   monotonic transformation of beta diversity)
#' }
#' @references
#' Szepannek, G. (2018) clustMixType: User-friendly clustering of mixed-type
#' data in R. \emph{The R Journal}, \strong{10}: 200-208.
#' @examples
#' \donttest{
#' library(terra)
#' data(africa)
#' p <- vect(system.file("ex/sa.json", package = "phyloregion"))
#' fb <- functional_beta(x=africa$comm, trait = africa$trait)
#' p <- phyloregion(fb[[1]], pol = p)
#' plot(p)
#' }
#' @export
functional_beta <- function(x, trait = NULL, bin = 10, na.rm = "no",
                            quick_elbow = FALSE, abundance = FALSE, ...) {
    if(is(x, "sparseMatrix")) {
        x <- sparse2long(x)
    }
    x$species <- gsub(" ", "_", x$species)

    trait <- .mtchtext(trait)
    trait <- as.data.frame(trait)
    trait <- trait[!duplicated(trait[ , "species"]),]
    trait$species <- gsub(" ", "_", trait$species)
    index <- intersect(x$species, trait$species)
    submat <- x[x$species %in% index, ]
    trait <- trait[trait$species %in% index, ]

    row.names(trait) <- trait$species
    trait <- trait[ , !colnames(trait) %in% "species"]

    x1 <- Filter(is.character, trait)
    x1[] <- lapply(x1, as.factor)

    x2 <- Filter(is.numeric, trait)
    x2[] <- apply(x2, 2, as.numeric)
    mm <- data.frame(x1, x2)

    # apply k-prototypes
    if(quick_elbow) {
        Es <- unlist(lapply(seq.int(2, bin), function(i) {
            kpres <- kproto(mm, k = i, verbose = FALSE, ...)
            y <- kpres$tot.withinss
            names(y) <- parent.frame()$i[]
            return(y)
        }))
        newk <- lbow(Es)$k
        mm <- kproto(mm, k = newk, na.rm = na.rm, ...)
    } else mm <- kproto(mm, k = bin, na.rm = na.rm, ...)

    memb <- mm$cluster

    foox <- function(x) {
        k <- max(x)
        res <- vector("list", k)
        for (i in seq_len(k)) {res[[i]] <- names(x)[x==i]}
        res
    }
    d <- foox(memb)
    m <- do.call(rbind, lapply(d, function(y) cbind(species=y,
                                                    cluster=parent.frame()$i[],
                                                    grouping="lumper",
                                                    ntax=length(y))))
    r <- data.frame(m[m[, "ntax"] > 1,])
    # MASTER_DISTINCT -- SPECIES WHICH ARE DISTINCT AND CAN'T BE LUMPED
    w <- submat[!submat$species %in% as.vector(r$species), ]
    spp <- unique(as.character(r$cluster))

    spo <- lapply(spp, function(y) {
        v <- r[r$cluster==y, ]
        v <- submat[submat$species %in% as.vector(v$species), ]
        cbind(v["grids"], species=v$species[1])
    })
    z <- rbind(w, do.call(rbind, spo))
    mx <- long2sparse(z)
    if(!abundance) mx@x[] <- 1
    res <- beta_diss(mx, ...)
    result <- list(res[[1]], res[[2]], res[[3]],
                   optimal_bins = if(quick_elbow) newk else bin)
    names(result)[1:3] <- c(names(res)[1], names(res)[2], names(res)[3])
    return(result)
}
