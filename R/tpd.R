# Functional alpha and beta diversity from trait probability densities (TPD)
# for mixed continuous and categorical traits.

# internally used in tpd_community
# mixed traits -> continuous axes (+ one group per combination of categories)
.tpd_axes <- function(trait, max.axes = 3, max.groups = 8, k = 2) {
  num <- Filter(is.numeric, trait)
  cat <- Filter(function(v) is.character(v) || is.factor(v), trait)
  ng <- if (ncol(cat)) nlevels(interaction(cat, drop = TRUE)) else 1L
  if (ncol(num) >= 1 && ncol(num) <= max.axes && ng <= max.groups) {
    axes <- scale(as.matrix(num))
    group <- if (ncol(cat)) interaction(cat, drop = TRUE, sep = "_") else
      factor(rep("all", nrow(trait)))
    method <- "direct"
  } else {
    Z <- if (ncol(num)) scale(as.matrix(num)) else NULL
    for (v in names(cat)) {
      D <- model.matrix(~ f - 1, data.frame(f = factor(cat[[v]])))
      p <- colMeans(D)
      Z <- cbind(Z, sweep(D, 2, p) / sqrt(1 - sum(p^2)))
    }
    pca <- prcomp(Z)
    axes <- pca$x[, seq_len(k), drop = FALSE]
    group <- factor(rep("all", nrow(trait)))
    method <- "pca"
  }
  rownames(axes) <- rownames(trait)
  list(axes = axes, group = group, method = method)
}


# internally used in tpd_alpha and tpd_beta_core
tpd_community <- function(x, trait, bins = 50, bw = NULL, quantile = 0.99,
                          cut = 3, ...) {
  ta <- .tpd_axes(trait, ...)
  a <- ta$axes
  d <- ncol(a)
  if (is.null(bw)) bw <- apply(a, 2, bw.nrd0)
  lo <- apply(a, 2, min) - cut * bw
  hi <- apply(a, 2, max) + cut * bw
  mids <- lapply(seq_len(d), function(j) seq(lo[j], hi[j], length.out = bins))
  step <- vapply(mids, function(m) m[2] - m[1], numeric(1))
  g <- as.integer(ta$group)
  K <- lapply(seq_len(nrow(a)), function(i) {
    z <- lapply(seq_len(d), function(j) {
      k <- which(abs(mids[[j]] - a[i, j]) <= cut * bw[j])
      list(k = k, w = dnorm(mids[[j]][k], a[i, j], bw[j]))
    })
    idx <- as.matrix(expand.grid(lapply(z, `[[`, "k")))
    w <- Reduce(`*`, expand.grid(lapply(z, `[[`, "w")))
    j <- as.vector(1 + (idx - 1) %*% bins^(0:(d - 1))) + (g[i] - 1) * bins^d
    w <- w / sum(w)
    # keep the smallest set of bins holding `quantile` of the density
    o <- order(w, decreasing = TRUE)
    o <- o[seq_len(which(cumsum(w[o]) >= quantile - 1e-12)[1])]
    list(j = j[o], x = w[o] / sum(w[o]))
  })
  K <- Matrix::sparseMatrix(as.integer(rep(seq_along(K), lengths(lapply(K, `[[`, "j")))),
                            as.integer(unlist(lapply(K, `[[`, "j"))),
                            x = unlist(lapply(K, `[[`, "x")),
                            dims = c(nrow(a), nlevels(ta$group) * bins^d),
                            dimnames = list(rownames(a), NULL))
  x <- x[, rownames(a)]
  x@x[] <- 1
  commtpd <- Diagonal(x = 1 / rowSums(x)) %*% x %*% K
  commtpd <- drop0(commtpd)
  rownames(commtpd) <- rownames(x)
  list(Matrix = commtpd, bin.volume = prod(step), K = K, method = ta$method)
}


# internally used in tpd_alpha and tpd_beta_core (like match_phylo_comm)
.match_trait_comm <- function(x, trait) {
  if (inherits(x, "matrix") && ncol(x) > 2) x <- Matrix(x, sparse = TRUE)
  if (!is(x, "sparseMatrix")) stop("x needs to be a sparse matrix!")
  trait <- .mtchtext(trait)
  trait$species <- gsub(" ", "_", trait$species)
  trait <- trait[!duplicated(trait$species), ]
  trait <- trait[complete.cases(trait), ]
  index <- intersect(colnames(x), trait$species)
  if (length(index) < 2) stop("Fewer than 2 species shared by x and trait!")
  rownames(trait) <- trait$species
  trait <- trait[index, setdiff(names(trait), "species"), drop = FALSE]
  x <- x[, index, drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  list(x = x, trait = trait)
}


#' Functional alpha diversity using trait probability densities
#'
#' \code{tpd_alpha} computes functional richness, evenness and redundancy of
#' each grid cell from the trait probability density (TPD) of its species,
#' for any mix of continuous and categorical traits.
#'
#' @param x an object of class Matrix or matrix, where rows are grid cells and
#' columns are species.
#' @param trait A data frame with the first column labeled \dQuote{species}
#' and the remaining columns containing continuous (numeric) and/or
#' categorical (character or factor) functional traits.
#' @param \dots Further arguments passed to the internal TPD builder:
#' \code{bins} (marks per trait axis, default 50), \code{bw} (kernel width
#' per standardized trait axis, default \code{bw.nrd0}), \code{quantile}
#' (share of each species' density kept, default 0.99, so that thin tails do
#' not count as occupied trait space), \code{max.axes} and \code{max.groups}
#' (with more continuous traits or more combinations of categories than
#' these, default 3 and 8, all traits are summarized by \code{k} = 2
#' mixed-PCA axes).
#' @details Continuous traits are standardized to mean 0 and standard
#' deviation 1, so transform skewed traits (e.g. log) beforehand. Each
#' combination of categorical traits gets its own trait space, so density
#' never spreads between categories.
#' @return A data frame with one row per grid cell: \code{SR} (species
#' richness), \code{FRic} (volume of trait space occupied), \code{FEve}
#' (evenness of the TPD), \code{FRed} (functional redundancy), community
#' weighted means of each continuous trait (\code{CWM_}) and the share of
#' each category (\code{share_}).
#' @references
#' Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2016) Traits without
#' borders: integrating functional diversity across scales. \emph{Trends in
#' Ecology & Evolution} \strong{31}: 382--394.
#'
#' Carmona, C.P., Tamme, R., Partel, M. et al. (2021) Erosion of global
#' functional diversity across the tree of life. \emph{Science Advances}
#' \strong{7}: eabf2675.
#' @seealso \code{\link{tpd_beta}}, \code{\link{PD}}
#' @examples
#' data(africa)
#' tr <- data.frame(species = africa$trait$species,
#'                  height = log10(africa$trait$h_max),
#'                  habit = africa$trait$habit)
#' fa <- tpd_alpha(africa$comm, tr)
#' head(fa)
#' @rdname tpd_alpha
#' @importFrom Matrix sparseMatrix Diagonal drop0 rowSums
#' @importFrom stats bw.nrd0 dnorm model.matrix prcomp complete.cases sd
#' @export
tpd_alpha <- function(x, trait, ...) {
  dat <- .match_trait_comm(x, trait)
  z <- tpd_community(dat$x, dat$trait, ...)
  T <- z$Matrix
  s <- as(T, "TsparseMatrix")
  s <- list(i = s@i + 1L, x = s@x)
  n_occ <- tabulate(s$i, nrow(T))
  pres <- dat$x
  pres@x[] <- 1
  Kb <- z$K
  Kb@x[] <- 1
  W <- Diagonal(x = 1 / rowSums(pres)) %*% pres
  res <- data.frame(grids = rownames(T), SR = rowSums(pres),
    FRic = n_occ * z$bin.volume,
    FEve = as.vector(tapply(pmin(s$x, 1 / n_occ[s$i]), s$i, sum)),
    FRed = pmax(0, as.vector(rowSums(T * (pres %*% Kb)) - 1)),
    row.names = NULL)
  num <- Filter(is.numeric, dat$trait)
  if (ncol(num)) {
    cwm <- as.matrix(W %*% as.matrix(num))
    colnames(cwm) <- paste0("CWM_", colnames(num))
    res <- cbind(res, cwm)
  }
  for (v in names(Filter(function(v) is.character(v) || is.factor(v), dat$trait))) {
    D <- model.matrix(~ f - 1, data.frame(f = factor(dat$trait[[v]])))
    sh <- as.matrix(W %*% D)
    colnames(sh) <- paste0("share_", sub("^f", "", colnames(D)))
    res <- cbind(res, sh)
  }
  res
}


#' Functional beta diversity using trait probability densities
#'
#' \code{tpd_beta_core} computes efficiently for large community matrices
#' the overlap between the trait probability densities (TPD) of pairs of grid
#' cells, and splits the non-overlapping part into the part in trait space
#' used by only one grid cell (unique) and the part in trait space used by
#' both (nested).
#'
#' @aliases tpd_beta tpd_beta_core
#' @param x an object of class Matrix or matrix, where rows are grid cells and
#' columns are species.
#' @param trait A data frame with the first column labeled \dQuote{species}
#' and the remaining columns containing continuous and/or categorical traits.
#' @param \dots Further arguments passed to the internal TPD builder, see
#' \code{\link{tpd_alpha}}.
#' @return \code{tpd_beta_core} returns an object of class "tpd.betapart",
#' a list of \code{dist} objects: \code{shared} (overlap), \code{not.shared}
#' (1 - overlap), \code{unique.i}, \code{unique.j}, \code{nested.i} and
#' \code{nested.j}.
#'
#' \code{tpd_beta} returns a list with three dissimilarity matrices:
#' \itemize{
#'   \item \code{beta.unique} dissimilarity due to trait space occupied by
#'   only one of the grid cells (turnover-like)
#'   \item \code{beta.nested} dissimilarity due to different densities in
#'   trait space shared by both grid cells (nestedness-like)
#'   \item \code{beta.total} total dissimilarity, 1 - overlap
#' }
#' @references
#' Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2016) Traits without
#' borders: integrating functional diversity across scales. \emph{Trends in
#' Ecology & Evolution} \strong{31}: 382--394.
#' @seealso \code{\link{tpd_alpha}}, \code{\link{phylobeta}},
#' \code{\link{beta_diss}}
#' @examples
#' data(africa)
#' tr <- data.frame(species = africa$trait$species,
#'                  height = log10(africa$trait$h_max),
#'                  habit = africa$trait$habit)
#' tbc <- tpd_beta_core(africa$comm, tr)
#' tb <- tpd_beta(africa$comm, tr)
#' @rdname tpd_beta
#' @importFrom Matrix tcrossprod colSums
#' @importFrom stats dist
#' @export
tpd_beta_core <- function(x, trait, ...) {
  dat <- .match_trait_comm(x, trait)
  z <- tpd_community(dat$x, dat$trait, ...)
  T <- z$Matrix[, colSums(z$Matrix) > 0]
  Labels <- rownames(T)
  l <- nrow(T)

  NOT.SHARED <- dist(as.matrix(T), method = "manhattan") / 2
  Tb <- T
  Tb@x[] <- 1
  OCC <- as.matrix(tcrossprod(T, Tb))   # density of row cell where column cell has some

  unique.j <- 1 - as.dist(t(OCC))
  unique.i <- 1 - as.dist(OCC)
  nested.j <- NOT.SHARED - unique.j
  nested.i <- NOT.SHARED - unique.i

  at <- list(Labels = Labels, Size = l, class = "dist", Diag = FALSE,
        Upper = FALSE)
  attributes(NOT.SHARED) <- at
  attributes(unique.i) <- at
  attributes(unique.j) <- at
  attributes(nested.i) <- at
  attributes(nested.j) <- at
  SHARED <- 1 - NOT.SHARED
  res <- list(shared = SHARED, not.shared = NOT.SHARED,
              unique.i = unique.i, unique.j = unique.j,
              nested.i = nested.i, nested.j = nested.j)
  class(res) <- "tpd.betapart"
  res
}


#' @rdname tpd_beta
#' @export
tpd_beta <- function(x, trait, ...) {
  res <- tpd_beta_core(x, trait, ...)
  mU <- pmin(res$unique.i, res$unique.j)
  PU <- 2 * mU / (res$nested.i + res$nested.j + 2 * mU)
  PU[res$not.shared == 0] <- 0
  at <- attributes(res$not.shared)
  beta.unique <- res$not.shared * PU
  beta.nested <- res$not.shared - beta.unique
  attributes(beta.unique) <- at
  attributes(beta.nested) <- at
  list(beta.unique = beta.unique, beta.nested = beta.nested,
       beta.total = res$not.shared)
}


#' Global trait spectra of groups of species
#'
#' \code{tpd_groups} turns a grouping of species (for example threatened
#' and non-threatened species, or plants of different habits) into a
#' groups x species matrix. Use it in place of a community matrix in
#' \code{\link{tpd_alpha}} and \code{\link{tpd_beta}} to describe the trait
#' spectrum of each group and how much the spectra overlap, as in Carmona et
#' al. (2021).
#'
#' @param group A vector of group labels, one per species, named by species.
#' @param all Logical, whether to add a row with all species together.
#' @return A sparse matrix with groups in rows and species in columns.
#' @references
#' Carmona, C.P., Tamme, R., Partel, M. et al. (2021) Erosion of global
#' functional diversity across the tree of life. \emph{Science Advances}
#' \strong{7}: eabf2675.
#' @seealso \code{\link{tpd_alpha}}, \code{\link{tpd_beta}},
#' \code{\link{tpd_extinction}}
#' @examples
#' data(africa)
#' tr <- data.frame(species = africa$trait$species,
#'                  height = log10(africa$trait$h_max),
#'                  habit = africa$trait$habit)
#' status <- ifelse(africa$IUCN$IUCN %in% c("CR", "EN", "VU"),
#'                  "threatened", "not_threatened")
#' names(status) <- africa$IUCN$Species
#' g <- tpd_groups(status)
#' tpd_alpha(g, tr)          # trait spectrum of each group
#' tpd_beta(g, tr)           # how different the spectra are
#' @rdname tpd_groups
#' @export
tpd_groups <- function(group, all = TRUE) {
  if (is.null(names(group))) stop("group needs species names!")
  group <- group[!is.na(group)]
  names(group) <- gsub(" ", "_", names(group))
  f <- factor(group)
  x <- Matrix::sparseMatrix(as.integer(f), seq_along(f), x = 1,
                            dimnames = list(levels(f), names(group)))
  if (all) {
    x <- rbind(x, all = 1)
  }
  x
}


#' Simulated extinctions and the loss of functional diversity
#'
#' \code{tpd_extinction} removes species at random according to their
#' extinction probabilities (for example from their IUCN category), many
#' times, and compares the trait spectrum before and after. As a null model,
#' the same number of species is removed completely at random.
#'
#' @param trait A data frame with the first column labeled \dQuote{species}
#' and the remaining columns containing continuous and/or categorical traits.
#' @param prob A numeric vector of extinction probabilities (0 to 1), one per
#' species, named by species.
#' @param nsim Number of simulations for each scenario.
#' @param \dots Further arguments passed to the internal TPD builder, see
#' \code{\link{tpd_alpha}}.
#' @return A list with \code{sims}, a data frame with one row per simulation
#' (scenario, species lost, FRic, FEve and FRed relative to the current
#' spectrum, and overlap with the current spectrum), and \code{summary}, the
#' mean and 95\% range of each quantity for the extinction scenario and for
#' random loss of the same number of species.
#' @references
#' Carmona, C.P., Tamme, R., Partel, M. et al. (2021) Erosion of global
#' functional diversity across the tree of life. \emph{Science Advances}
#' \strong{7}: eabf2675.
#' @seealso \code{\link{tpd_groups}}, \code{\link{tpd_beta}}
#' @examples
#' data(africa)
#' tr <- data.frame(species = africa$trait$species,
#'                  height = log10(africa$trait$h_max),
#'                  habit = africa$trait$habit)
#' # hypothetical extinction probabilities for each IUCN category
#' p <- c(LC = 0.001, NT = 0.01, VU = 0.1, EN = 0.5, CR = 0.9)
#' prob <- setNames(p[africa$IUCN$IUCN], africa$IUCN$Species)
#' set.seed(1)
#' ex <- tpd_extinction(tr, prob, nsim = 20)
#' ex$summary
#' @rdname tpd_extinction
#' @importFrom stats runif
#' @export
tpd_extinction <- function(trait, prob, nsim = 100, ...) {
  names(prob) <- gsub(" ", "_", names(prob))
  prob <- prob[!is.na(prob)]
  tm <- .mtchtext(trait)
  tm <- tm[complete.cases(tm), ]
  sp <- intersect(names(prob), gsub(" ", "_", tm$species))
  if (length(sp) < 3) stop("Fewer than 3 species shared by prob and trait!")
  prob <- prob[sp]
  S <- length(sp)
  keep_ext <- lapply(seq_len(nsim), function(i) which(runif(S) > prob))
  keep_rnd <- lapply(keep_ext, function(k) sample(S, length(k)))
  keep <- c(list(seq_len(S)), keep_ext, keep_rnd)
  lab <- c("current", paste0("extinction_", seq_len(nsim)),
           paste0("random_", seq_len(nsim)))
  x <- Matrix::sparseMatrix(rep(seq_along(keep), lengths(keep)),
                            unlist(keep), x = 1,
                            dims = c(length(keep), S),
                            dimnames = list(lab, sp))
  fa <- tpd_alpha(x, trait, ...)
  tb <- tpd_beta(x, trait, ...)
  cur <- fa[fa$grids == "current", ]
  sims <- fa[fa$grids != "current", ]
  sims <- data.frame(sim = sims$grids,
    scenario = sub("_[0-9]+$", "", sims$grids),
    lost = S - sims$SR,
    FRic = sims$FRic / cur$FRic,
    FEve = sims$FEve / cur$FEve,
    FRed = sims$FRed / cur$FRed,
    overlap = 1 - as.matrix(tb$beta.total)["current", sims$grids],
    row.names = NULL)
  vars <- c("lost", "FRic", "FEve", "FRed", "overlap")
  summ <- do.call(rbind, lapply(split(sims, sims$scenario), function(d) {
    data.frame(scenario = d$scenario[1], quantity = vars,
      mean = vapply(vars, function(v) mean(d[[v]]), numeric(1)),
      lower = vapply(vars, function(v) unname(quantile(d[[v]], 0.025)), numeric(1)),
      upper = vapply(vars, function(v) unname(quantile(d[[v]], 0.975)), numeric(1)),
      row.names = NULL)
  }))
  rownames(summ) <- NULL
  list(sims = sims, summary = summ)
}
