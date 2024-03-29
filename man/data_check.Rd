% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_check.R
\name{data_check}
\alias{data_check}
\title{Anomaly checking for network node variables}
\usage{
data_check(covariates, groups)
}
\arguments{
\item{covariates}{A numeric sample matrix with n rows (all samples) and p columns of network covariates}

\item{groups}{A character vector indicating a group in which each sample belongs to.}
}
\value{
The first table provides information about the empirical distributions of each variable:
Min, max, skewness, and kurtosis. The second table shows group-specific sample sizes.
}
\description{
Anomaly checking for network node variables
}
\examples{
print("anomaly_check(covariates)")
}
