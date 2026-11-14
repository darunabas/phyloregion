## Submission

This is an update from version 1.0.9 to 2.1.0.

* Fixed `sdm()` failing with `algorithm = "all"` (the default) and `algorithm = "GLM"`.
* Replaced the deprecated `structure(.Names = )` call that produced a NOTE on r-devel.
* New functions for functional diversity based on trait probability densities
  (`tpd_alpha()`, `tpd_beta_core()`, `tpd_beta()`, `tpd_groups()`, `tpd_extinction()`)
  and a new vignette.

## Test environments

* local macOS (aarch64), R 4.6.0
* win-builder: R-devel and R-release (Windows)

## R CMD check results

0 errors | 0 warnings | 0 notes on local macOS.

On win-builder there is 1 NOTE: "Possibly misspelled words in DESCRIPTION: Carmona". This is an author name in the cited reference (Carmona et al. 2016) and is spelled correctly.

## Reverse dependencies

There is one reverse dependency on CRAN, canaper (1.0.1). It was checked with R CMD check against this version: Status OK, no new problems.
