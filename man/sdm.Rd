% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdm.R
\name{sdm}
\alias{sdm}
\title{Species distribution models for a range of algorithms}
\usage{
sdm(
  x,
  pol = NULL,
  predictors = NULL,
  blank = NULL,
  res = 1,
  tc = 2,
  lr = 0.001,
  bf = 0.75,
  n.trees = 50,
  step.size = n.trees
)
}
\arguments{
\item{x}{A dataframe containing the species occurrences
and geographic coordinates. Column 1 labeled as "species", column 2 "lon",
column 3 "lat".}

\item{pol}{A polygon shapefile specifying the boundary to restrict the
prediction. If not specified, a minimum convex polygon is estimated using
the input data frame of species occurrences.}

\item{predictors}{RasterStack of environmental descriptors on which
the models will be projected}

\item{blank}{A blank raster upon which the prediction layer is aggregated to.}

\item{res}{Desired resolution of the predicted potential species distribution
(if blank raster is not specified).}

\item{tc}{Integer. Tree complexity. Sets the complexity of individual trees}

\item{lr}{Learning rate. Sets the weight applied to individual trees}

\item{bf}{Bag fraction. Sets the proportion of observations used in selecting variables}

\item{n.trees}{Number of initial trees to fit. Set at 50 by default}

\item{step.size}{Number of trees to add at each cycle}
}
\value{
A list with the following objects:
\itemize{
  \item \code{ensemble_raster} The ensembled raster that predicts
  the potential species distribution.
  \item \code{ensemble_AUC} The median AUCs of models.
  \item \code{data} The dataframe that was used to implement the model.
  \item \code{indiv_models} Raster layers for the separate models that
  predict the potential species distribution.
  \item \code{single_AUCs} The AUCs for the seperate models.
}
}
\description{
This function computes species distribution models using
four modelling algorithms: generalized linear models,
generalized boosted models, random forests, and maximum entropy.
}
\references{
Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006) Maximum entropy
modeling of species geographic distributions. \emph{Ecological Modelling}
\strong{190}: 231-259.
}