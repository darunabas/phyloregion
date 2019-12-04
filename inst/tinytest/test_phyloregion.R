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

pd_bioregion <- PD.sparse(com, tree)

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
