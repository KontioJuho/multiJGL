% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiJGL.R
\name{multiJGL}
\alias{multiJGL}
\title{Linear and nonlinear multiclass network estimation with joint regularization over categorical groups}
\usage{
multiJGL(
  node.covariates = node.covariates,
  grouping.factor = grouping.factor,
  penalty.lin = "fused",
  penalty.nonlin = "fused",
  lin_lambda1 = 0.1,
  lin_lambda2 = 0.025,
  nonlin_lambda1 = 0.1,
  nonlin_lambda2 = 0.025,
  tol.linear = 1e-05,
  tol.nonlinear = 1e-05,
  subsample.pseudo_obs = FALSE,
  omit.rate = 2L,
  ...
)
}
\arguments{
\item{node.covariates}{An nxp dimensional matrix of p covariates measured over n samples.}

\item{grouping.factor}{A grouping factor for creating observational classes.}

\item{penalty.lin}{Specify "fused" or "group" penalty type for the linear JGL algorithm.}

\item{penalty.nonlin}{Specify "fused" or "group" penalty type for the nonlinear JGL algorithm.
#See the explanation for the fused and group penalties in the JGL package
#The original JGL CRAN repository: https://CRAN.R-project.org/package=JGL
#The following penalty parameters are given in pairs to --
#separately assign the amount of regularizations for linear and nonlinear parts}

\item{lin_lambda1}{The l1-penalty parameter for the linear JGL to regulate within group network densities}

\item{lin_lambda2}{The l1-penalty parameter for the nonlinear JGL.}

\item{nonlin_lambda1}{The fusion penalty parameter for the linear JGL.}

\item{nonlin_lambda2}{The fusion penalty parameter for the nonlinear JGL.}

\item{tol.linear}{Convergence criterion for the linear part (see the JGL package for details).}

\item{tol.nonlinear}{Convergence criterion for the nonlinear part.
#Subsampling procedure over pseudo-observations if the number of observations is large already in the original sets.}

\item{subsample.pseudo_obs}{Should the subsampling procedure be used over the pseudo-observations.}

\item{omit.rate}{An integer: Omit rate for the subsampling pcocedure between 2L and 5L}

\item{...}{Additional parameter for the nonlinear JGL.}
}
\description{
Linear and nonlinear multiclass network estimation with joint regularization over categorical groups
}
\examples{
 print("net <- multiJGL(node.covariates, grouping.factor)")
}
