# bioregion
This R package is for mapping biodiversity and for conservation. It is simple, fast, and particularly tailored for handling large datasets.
### Authors
[Barnabas Daru](https://barnabasdaru.com/) 

[Klaus Schliep](https://kschliep.netlify.com/)
### Citation
If you find ```bioregion``` helpful, please cite as:

Daru B. H. & Schliep, K. bioregion: an R package for mapping biodiversity and conservation (R package version 0.1.0. https://github.com/darunabas/bioregion, 2019).
# Introduction
This tutorial is an introduction to using `R` for analysing geographic data in biodiversity science and conservation. Specifically, you will be testing my new `R` package called `bioregion` which is still in a beta-testing phase. The `bioregion` package is a tool for mapping various facets of biodiversity ranging from local (alpha-) to between community (beta-) diversity.

The `bioregion` package will introduce the basics of mapping various facets of spatial data ranging from species richness, endemism, to threat as evaluated by the International Union for the Conservation of Nature as well as beta diversity metrics. More advanced implementations of `bioregion` is the addition of phylogenetic information to quantify evolutionary diversity including phylogenetic diversity, phylogenetic endemism, and evolutionary distinctiveness and global endangerment.

A major feature of `bioregion` is its ability to handle large datasets spanning 1000s to hundreds of thousands of taxa and spanning large geographic extents.
# Installation
The `bioregion` package is available from github. First, you will need to install the `devtools` package. In R, type:
```{r, echo=TRUE}
#install.packages("devtools") # uncommenting this will install the package
```
Next, load the `devtools` package.
```{r, message=FALSE, results='hide', warning=FALSE}
library(devtools)
```
Then install the `bioregion` package from github:
```{r, echo=TRUE}
#install_github("darunabas/bioregion") # uncommenting this will install bioregion package
```
Load the `bioregion` package:
```{r, echo=TRUE}
library(bioregion)
```
Although the package's strong focus is for mapping biodiversity patterns, we will draw from other packages including: `raster`, `Matrix`, `ape`, `data.table` and `rgeos`.

```{r, message=FALSE, results='hide', warning=FALSE}
z <- c("raster", "Matrix", "ape", "colorRamps", "data.table", "rgeos")
# install.packages(z) # uncommenting this will install the packages
lapply(z, library, character.only = TRUE) # load the required packages
```