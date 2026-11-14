# [November 13, 2026] Changes in `phyloregion` version 2.1.0.

* Revised and improved the function sdm to return the evaluation threshold.

* Fixed sdm() failing with algorithm = "all" (the default) and algorithm = "GLM".

* sdm() now returns the continuous habitat suitability for every algorithm in raw_habitatsuit.

* Replaced the deprecated structure(.Names = ) call in phylobeta_core (CRAN check NOTE).

* New functions tpd_alpha(), tpd_beta_core() and tpd_beta() compute functional alpha diversity (richness, evenness, redundancy, community weighted means) and functional beta diversity split into unique (turnover-like) and nested components, using trait probability densities for any mix of continuous and categorical traits (Carmona et al. 2016).

* New functions tpd_groups() and tpd_extinction() describe the trait spectra of groups of species (for example threatened versus non-threatened) and simulate the erosion of functional diversity under extinction scenarios, with random loss of the same number of species as a null model (Carmona et al. 2021).

* New vignette: Functional diversity with trait probability densities.

# [July 14, 2025] Changes in `phyloregion` version 1.0.9.

* Prepared and updated package for CRAN.

* Revised and updated plot_pie to plot_spatial_membership and now uses blended color visualization.

* Automatically adjusts radius depending on the projection system (e.g., EPSG:4326 vs EPSG:6933).

* Added new functions backg, randpoints, rast_quantile, rast_hotspot

* Updated documentation and examples accordingly.


# [April 25, 2023] Changes in `phyloregion` version 1.0.8.

* Prepared and updated package for CRAN.

* Revised all spatial functions to use functions in `terra` throughout.

* Discontinued use of `raster`, `rgeos` and `sp` packages.

* pd and beta diversity functions now accepts a community matrix of class matrix
  
  or an `phyloseq` object as input.


# [April 30, 2021] Changes in `phyloregion` version 1.0.5

* Revised `plot_swatch` to use `col` instead of palette argument.

* Added new function, `sdm`, for species distribution modelling.

* Improved functions for transforming converting input data.


# [July 19, 2020] Changes in `phyloregion` version 1.0.3

* Removed `plot_evoldistinct` because the functionality is replaced 
with `plot_swatch`.

* `plot.phyloregion` and the generic `plot` function are now exported and both 
produce the same result. 

* Changed title to "`phyloregion`: R package for biogeographic 
regionalization and macroecology"

* Improvements in the vignettes with several examples.


# [March 30, 2020] `phyloregion` version 1.0.0

* This is the first CRAN release of phyloregion.
