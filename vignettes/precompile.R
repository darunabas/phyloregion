# Precompiled vignettes that depend on API key
# Must manually move images from phyloregion/ to phyloregion/vignettes/ after knit
#

library(knitr)
knit("vignettes/Benchmark.Rmd.orig", "vignettes/Benchmark.Rmd")

knit("vignettes/scalability.Rmd.orig", "vignettes/scalability.Rmd")

# run pkgdown::build_site() afterwards
