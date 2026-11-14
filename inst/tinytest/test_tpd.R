## TPD functional diversity: checks against the hand example
## 4 species, 4 grid cells; height (continuous) + growth form (categorical)

library(phyloregion)
library(Matrix)

tr <- data.frame(species = c("A", "B", "C", "D"),
                 logh = log10(c(18, 26, 0.4, 0.7)),
                 growth_form = c("tree", "tree", "herb", "herb"))
comm <- Matrix(c(1,1,0,0,    # v1: A B     (two trees)
                 1,0,1,0,    # v2: A C     (a tree and a herb)
                 1,1,1,1,    # v3: A B C D (all four)
                 0,0,1,1),   # v4: C D     (two herbs)
               nrow = 4, byrow = TRUE, sparse = TRUE,
               dimnames = list(paste0("v", 1:4), c("A", "B", "C", "D")))
bw <- 0.15 / sd(tr$logh)     # bell width of 0.15 log10 m, as in the hand example

# quantile = 1 reproduces the hand example, which cut each bell at 3 widths only
tb <- suppressMessages(tpd_beta(comm, tr, bw = bw, bins = 400, quantile = 1))
fa <- suppressMessages(tpd_alpha(comm, tr, bw = bw, bins = 400, quantile = 1))
get <- function(d, a, b) as.matrix(d)[a, b]

## beta: structure
expect_true(all(sapply(tb, inherits, "dist")))
expect_equal(labels(tb$beta.total), paste0("v", 1:4))
expect_equal(as.vector(tb$beta.unique + tb$beta.nested), as.vector(tb$beta.total))
expect_true(all(tb$beta.total >= 0 & tb$beta.total <= 1 + 1e-12))

## beta: values worked out by hand
# v1 (trees) vs v4 (herbs): no shared trait space -> all unique (turnover)
expect_equal(get(tb$beta.total, "v1", "v4"), 1)
expect_equal(get(tb$beta.unique, "v1", "v4"), 1)
# v1 and v2 share exactly half a bucket (species A)
expect_equal(get(tb$beta.total, "v1", "v2"), 0.5)
# v1 sits inside v3 (v3 = v1 + herbs) -> all nested
expect_equal(get(tb$beta.total, "v1", "v3"), 0.5)
expect_equal(get(tb$beta.unique, "v1", "v3"), 0)
expect_equal(get(tb$beta.unique, "v3", "v4"), 0)
# v1 vs v2: mostly nested, a little unique (the upper tail of B)
expect_true(get(tb$beta.unique, "v1", "v2") < 0.05)
expect_equal(get(tb$beta.nested, "v1", "v2"), 0.48, tolerance = 0.02)
expect_equal(get(tb$beta.total, "v2", "v3"), 0.25, tolerance = 0.02)

## beta: finer ruler marks give the same answer (height is continuous)
tb2 <- suppressMessages(tpd_beta(comm, tr, bw = bw, bins = 800, quantile = 1))
expect_equal(as.vector(tb2$beta.total), as.vector(tb$beta.total), tolerance = 0.01)

## default 99% cutoff: the parts that do not depend on the tails stay the same
tq <- suppressMessages(tpd_beta(comm, tr, bw = bw, bins = 400))
expect_equal(get(tq$beta.total, "v1", "v4"), 1)
expect_equal(get(tq$beta.unique, "v1", "v4"), 1)
expect_equal(get(tq$beta.total, "v1", "v2"), 0.5)
expect_equal(get(tq$beta.unique, "v1", "v3"), 0)
expect_equal(as.vector(tq$beta.unique + tq$beta.nested), as.vector(tq$beta.total))

## alpha: values worked out by hand
expect_equal(fa$grids, paste0("v", 1:4))
expect_equal(fa$SR, c(2, 2, 4, 2))
expect_equal(fa$share_tree, c(1, 0.5, 0.5, 0))
expect_equal(round(10^fa$CWM_logh, 1), c(21.6, 2.7, 3.4, 0.5))
expect_equal(fa$FRed[2], 0)                 # tree and herb never share trait space
expect_true(all(fa$FRed[c(1, 3, 4)] > 0.8))
expect_true(fa$FRic[2] > fa$FRic[1])        # tree + herb covers more than two trees
expect_true(all(fa$FEve > 0 & fa$FEve <= 1))

## inputs: dense matrix and spaces in species names work too
tr_sp <- tr
tr_sp$species <- c("sp A", "sp B", "sp C", "sp D")
comm_sp <- as.matrix(comm)
colnames(comm_sp) <- c("sp_A", "sp_B", "sp_C", "sp_D")
tb3 <- suppressMessages(tpd_beta(comm_sp, tr_sp, bw = bw, bins = 400, quantile = 1))
expect_equal(as.vector(tb3$beta.total), as.vector(tb$beta.total))
expect_error(suppressMessages(tpd_beta(data.frame(a = 1), tr)))

## many traits: switches to mixed PCA and still returns valid densities
set.seed(1)
trm <- data.frame(species = c("A", "B", "C", "D"), t1 = rnorm(4), t2 = rnorm(4),
                  t3 = rnorm(4), t4 = rnorm(4), g = c("x", "y", "x", "y"))
tbm <- suppressMessages(tpd_beta(comm, trm))
expect_true(all(tbm$beta.total >= 0 & tbm$beta.total <= 1 + 1e-12))

## global spectra of groups (threatened vs not) -----------------------------
grp <- c(A = "tree", B = "tree", C = "herb", D = "herb")
g <- tpd_groups(grp)
expect_equal(rownames(g), c("herb", "tree", "all"))
expect_equal(as.vector(rowSums(g)), c(2, 2, 4))
tg <- suppressMessages(tpd_beta(g, tr, bw = bw, bins = 400, quantile = 1))
expect_equal(get(tg$beta.total, "herb", "tree"), 1)      # trees and herbs never overlap
expect_equal(get(tg$beta.total, "tree", "all"), 0.5)     # trees are half of all species
expect_equal(get(tg$beta.unique, "tree", "all"), 0)      # ... and sit inside it (nested)
expect_error(tpd_groups(c("a", "b")))                     # species names are required

## extinction scenarios -------------------------------------------------------
set.seed(2)
ex <- suppressMessages(tpd_extinction(tr, c(A = 0, B = 0, C = 1, D = 0), nsim = 5,
                                      bw = bw, bins = 400))
expect_true(all(ex$sims$lost == 1))                       # C always goes extinct
expect_equal(ex$sims$lost[ex$sims$scenario == "random"],
             ex$sims$lost[ex$sims$scenario == "extinction"])
expect_true(all(ex$sims$overlap >= 0 & ex$sims$overlap <= 1 + 1e-12))
expect_equal(sort(unique(ex$summary$scenario)), c("extinction", "random"))
