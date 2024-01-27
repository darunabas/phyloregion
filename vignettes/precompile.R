# Precompiled vignettes that depend on API key
# Must manually move images from phyloregion/ to phyloregion/vignettes/ after knit
#

library(knitr)
knit("vignettes/pd.Rmd.orig", "vignettes/pd.Rmd")

knit("vignettes/scalability.Rmd.orig", "vignettes/scalability.Rmd")

knit("vignettes/heterogeneity.Rmd.orig", "vignettes/heterogeneity.Rmd")

knit("vignettes/ed.Rmd.orig", "vignettes/ed.Rmd")

# run pkgdown::build_site() afterwards
# ch <- check_for_cran()
