# Precompiled vignettes that depend on API key
# Must manually move image files from phyloregion/ to phyloregion/vignettes/ after knit

library(knitr)
knit("vignettes/Benchmark.Rmd.orig", "vignettes/Benchmark.Rmd")
