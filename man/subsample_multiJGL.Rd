% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsample_multiJGL.R
\name{subsample_multiJGL}
\alias{subsample_multiJGL}
\title{Function for running multiJGL algorithm across subsamples and different sparsity levels}
\usage{
subsample_multiJGL(
  node.covariates = node.covariates,
  grouping.factor = grouping.factor,
  lambda.seq,
  by_value,
  num_repetitions = 10,
  penalty.lin = "fused",
  penalty.nonlin = "fused",
  lin_lambda1 = lambda.seq1[j],
  lin_lambda2 = 0.025,
  nonlin_lambda1 = lambda.seq1[j],
  nonlin_lambda2 = 0.025,
  tol.linear = 1e-05,
  tol.nonlinear = 1e-05,
  subsample.pseudo_obs = FALSE,
  omit.rate = 2L
)
}
\arguments{
\item{node.covariates}{An nxp dimensional matrix of p covariates measured over n samples.}

\item{grouping.factor}{A grouping factor for creating observational classes.}

\item{lambda.seq}{The range of different lambda values}

\item{by_value}{Specifies the density of lambda grid}

\item{num_repetitions}{The number of subsample analyses}

\item{penalty.lin}{Specify "fused" or "group" penalty type for the linear JGL algorithm.}

\item{penalty.nonlin}{Specify "fused" or "group" penalty type for the nonlinear JGL algorithm.}

\item{lin_lambda1}{The l1-penalty parameter for the linear JGL to regulate within group network densities}

\item{lin_lambda2}{The l1-penalty parameter for the nonlinear JGL.}

\item{nonlin_lambda1}{The fusion penalty parameter for the linear JGL.}

\item{nonlin_lambda2}{The fusion penalty parameter for the nonlinear JGL.}

\item{tol.linear}{Convergence criterion for the linear part (see the JGL package for details).}

\item{tol.nonlinear}{Convergence criterion for the nonlinear part.}

\item{subsample.pseudo_obs}{Should the subsampling procedure be used over the pseudo-observations.}

\item{omit.rate}{An integer: Omit rate for the subsampling pcocedure between 2L and 5L}
}
\description{
Function for running multiJGL algorithm across subsamples and different sparsity levels
}
\examples{
print("")
}
