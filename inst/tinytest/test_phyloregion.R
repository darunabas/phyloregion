## generate data

library(phyloregion)
library(betapart)
library(ape)
library(Matrix)

tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
com <- matrix(c(1,0,1,1,0,0,
                1,0,0,1,1,0,
                1,1,1,1,1,1,
                0,0,1,1,0,1), 6, 4,
              dimnames=list(paste0("g",1:6), tree$tip.label))

com_sparse <- Matrix(com, sparse=TRUE)

pd_bioregion <- PD(com_sparse, tree)

pbc_phyloregion <- phylobeta_core(com_sparse, tree)
pbc_betapart <- phylo.betapart.core(com, tree)

# check if betapart and phyloregion are equal
expect_equal(phylo.beta.multi(pbc_phyloregion),
             phylo.beta.multi(pbc_betapart))
expect_equivalent(phylo.beta.pair(pbc_phyloregion),
                  phylo.beta.pair(pbc_betapart))

# test pd compare with picante
# test
pb_phyloregion <- phylobeta(com_sparse, tree)



data(africa)
long <- africa$comm

sparse <- sampl2sparse(long)
long2 <- sparse2sampl(sparse)

expect_equal(nrow(long2), nrow(long))

set.seed(42)
ind <- sample(60823, 1000, replace = TRUE)
long_ind <- long[ind,]

sparse_ind <- sampl2sparse(long_ind, method="nonphylo")
long_ind2 <- sparse2sampl(sparse_ind)

expect_equal(nrow(long_ind), nrow(long_ind2))


#check if beta_core is equal to betapart.core
M_sparse <- sampl2sparse(africa$comm, method = "nonphylo")
M_dense <- as(M_sparse, "matrix")

bc_dense <- betapart.core(M_dense)
bc_sparse <- beta_core(M_sparse)
# first element is data matrix, so we ignore it
expect_equal(bc_sparse[-1], bc_dense[-1])


# choropleth works (expect high positive correlation)
x <- rnorm(1000)
expect_true(cor(x, phyloregion:::choropleth3(x)) > 0.9)



# evol_distinct vs evol.distinct
if(requireNamespace("picante")){
  tree <- rcoal(100)
  expect_equivalent(evol_distinct(tree), picante::evol.distinct(tree)[,2])
  expect_equivalent(evol_distinct(tree, "fair"),
                    picante::evol.distinct(tree, "fair")[,2])
}

