## generate data

library(phyloregion)
library(betapart)
library(ape)

tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
com <- matrix(c(1,0,1,1,0,0,
                1,0,0,1,1,0,
                1,1,1,1,1,1,
                0,0,1,1,0,1), 6, 4,
              dimnames=list(paste0("g",1:6), tree$tip.label))

pd_bioregion <- PD(com, tree)

pbc_phyloregion <- phylobeta_core(com, tree)
pbc_betapart <- phylo.betapart.core(com, tree)

# check if betapart and phyloregion are equal
expect_equal(phylo.beta.multi(pbc_phyloregion),
             phylo.beta.multi(pbc_betapart))
expect_equivalent(phylo.beta.pair(pbc_phyloregion),
                  phylo.beta.pair(pbc_betapart))

# test pd compare with picante
# test
pb_phyloregion <- phylobeta(com, tree)



data(africa)
long <- africa$comm

sparse <- sampl2sparse(long)
long2 <- sparse2sampl(sparse)

expect_equal(nrow(long2), nrow(long))

set.seed(42)
ind <- sample(60823, 1000, replace = TRUE)
long_ind <- long[ind,]
#sparse_ind <- sampl2sparse(long_ind)
sparse_ind <- sampl2sparse(long_ind, method="nonphylo")
long_ind2 <- sparse2sampl(sparse_ind)

expect_equal(nrow(long_ind), nrow(long_ind2))




