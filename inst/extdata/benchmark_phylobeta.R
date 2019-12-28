rm(list = ls())

x <- c(100, 200, 500, 1000, 2000, 3000, 5000, 10000)

m <- matrix(0, length(x), 3,
                  dimnames=list(as.character(x),
                                c("picante", "betapart","phyloregion")))

if(require(ape) && require(phyloregion) && require(betapart) && require(picante)){
    i <- 100
    for (i in x)
    {
        # random tree
        tree <- as(rtree(i), "phylo")

        # random community matrix:
        com <- diag(1, i, i)
        colnames(com) <- tree$tip.label
        picante.time <- system.time(picante <- picante::phylosor(com, tree))["elapsed"]

        betapt.time <- system.time(betapart <- phylo.betapart.core(com, tree))["elapsed"]

        M <- comm2sparse(com)
        phyloregion.time <- system.time(phyloregion <- phylobeta(M, tree))["elapsed"]

        m[as.character(i),] <- c(picante.time, betapt.time, phyloregion.time)
    }
}

print(m)


C <- hcl.colors(4, "Zissou 1")

postscript("/Users/barnabasdaru/Dropbox/Projects/phyloregion_MS/manuscript/2019/12.December/Dec23/benchmark_phylobeta.eps", height = 8, width = 8)
plot(x, m[,"picante"], pch=21, bg=C[1],
     xlab="Phylogenetic tree size (number of tips)",
     ylab="Execution time (seconds)",
     ylim=c(0, max(m)), las=1, cex=2)
lines(x, m[,"picante"], lwd=2, col=C[1])

lines(x, m[,"betapart"], lwd=2, col=C[2])
points(x, m[,"betapart"], pch=21, bg=C[2], cex=2)

lines(x, m[,"phyloregion"], lwd=2, col=C[4])
points(x, m[,"phyloregion"], pch=21, bg=C[4], cex=2)


legend("topleft", c("picante", "betapart", "phyloregion"),
       pch=c(21, 21, 21), pt.bg=c(C, pt.cex=2, bty="n"), bty = "n")

dev.off()



cat("*** End of benchmarker.\n")
