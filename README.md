[![Build Status](https://travis-ci.org/darunabas/phyloregion.svg?branch=master)](https://travis-ci.org/darunabas/phyloregion)
[![codecov](https://codecov.io/gh/darunabas/phyloregion/branch/master/graph/badge.svg)](https://codecov.io/gh/darunabas/phyloregion)

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

# phyloregion <img src='man/figures/logo.png' align="right" width="120" />

This is a computational infrastructure for biogeographic regionalization (the classification of geographical areas in terms of their biotas) and spatial conservation in the R scientific computing environment. Previously it was only possible to perform analysis of biogeographic regionalization on small datasets, often using tools that are difficult to replicate. With macroecological datasets of ever increasing size and complexity, `phyloregion` offers the possibility of handling and executing large scale biogeographic regionalization efficiently and with extreme speed. It also allows fast and efficient analysis of more standard conservation measures such as phylogenetic diversity, phylogenetic endemism, evolutionary distinctiveness and global endangerment. `phyloregion` can run on any operating system (Mac, Linux, Windows or even high performance computing cluster) with `R` 3.6.0 (or higher) installed.
### Authors
[Barnabas Daru](https://barnabasdaru.com/) 

[Klaus Schliep](https://kschliep.netlify.com/)
### How to cite
The original implementation of ```phyloregion``` is described in:

Daru B.H., Karunarathne, P. & Schliep, K. (2020) phyloregion: R package for biogeographic regionalization and spatial conservation. **_bioRxiv_** 2020.02.12.945691 doi: 10.1101/2020.02.12.945691

It is based on the method described in:

Daru, B.H., Farooq, H., Antonelli, A. & Faurby, S. (2020) Endemism patterns are scale dependent. **_Coming soon_**. 

The original conceptual is described in:

Daru, B. H., Elliott, T. L., Park, D. S. & Davies, T. J. (2017), Understanding the processes underpinning patterns of phylogenetic regionalization. **_Trends in Ecology and Evolution_** __32__: 845-860. doi: 10.1016/j.tree.2017.08.013

### Feedback
If you have any questions, suggestions or issues regarding the package, please add them to [GitHub issues](https://github.com/darunabas/phyloregion/issues)

# Installation

`phyloregion` is an open-source and free package hosted on GitHub. You will need to install the `devtools` package. In R, type:

```
if (!requireNamespace("devtools", quietly = TRUE)) 
    install.packages("devtools") 
```

Then:

```
devtools::install_github("darunabas/phyloregion")
```

Load the `phyloregion` package:

```
library(phyloregion)
```
