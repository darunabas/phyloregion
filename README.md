# phyloregion
This R package is for analysis of biogeographic regionalization (the classification of geographical areas in terms of their biotas) and spatial conservation in the R scientific computing environment. Previously it was only possible to perform analysis of biogeographic regionalization on small datasets, often using tools that are difficult to replicate. With macroecological datasets of ever increasing size and complexity, `phyloregion` now offers the possibility of clustering large-scale species distributions combined with phylogenetic information using sparse matrix, determining optimal number of clusters, quantifying evolutionary distinctiveness of bioregions, and visualizing patterns and underlying processes. It is simple, fast, and particularly tailored for handling large datasets.
### Authors
[Barnabas Daru](https://barnabasdaru.com/) 

[Klaus Schliep](https://kschliep.netlify.com/)
### Citation
If you find ```phyloregion``` helpful, please cite as:

Daru B. H., Karunarathne, P. & Schliep, K. phyloregion: R package for biogeographic regionalization and spatial conservation (R package version 0.1.0. https://github.com/darunabas/phyloregion, 2019).
# Introduction
This tutorial is an introduction to using `R` for analyzing geographic data in biodiversity science and conservation. The `phyloregion` package is a tool for mapping various facets of biodiversity ranging from local (alpha-) to between community (beta-) diversity.

The `phyloregion` package will introduce the basics of mapping various facets of spatial data ranging from species richness, endemism, to threat as evaluated by the International Union for the Conservation of Nature as well as beta diversity metrics. More advanced implementations of `phyloregion` is the addition of phylogenetic information to quantify evolutionary diversity including phylogenetic diversity, phylogenetic endemism, and evolutionary distinctiveness and global endangerment.

A major feature of `phyloregion` is its ability to handle large datasets spanning 1000s to hundreds of thousands of taxa and spanning large geographic extents.
# Installation
The `phyloregion` package is available from github. First, you will need to install the `devtools` package. In R, type:
```{r, echo=TRUE}
#install.packages("devtools") # uncommenting this will install the package
```
Next, load the `devtools` package.
```{r, message=FALSE, results='hide', warning=FALSE}
library(devtools)
```
Then install the `phyloregion` package from github:
```{r, echo=TRUE}
#install_github("darunabas/phyloregion") # uncommenting this will install phyloregion package
```
Load the `phyloregion` package:
```{r, echo=TRUE}
library(phyloregion)
```
Although the package's strong focus is for mapping biodiversity patterns, we will draw from other packages including: `raster`, `Matrix`, `ape`, `data.table` and `rgeos`.

```{r, message=FALSE, results='hide', warning=FALSE}
z <- c("raster", "Matrix", "ape", "colorRamps", "data.table", "rgeos")
# install.packages(z) # uncommenting this will install the packages
lapply(z, library, character.only = TRUE) # load the required packages
```
