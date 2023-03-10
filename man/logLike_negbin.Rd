% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{logLike_negbin}
\alias{logLike_negbin}
\title{Log-likelihood for negative binomial with interfered responses}
\usage{
logLike_negbin(beta, alpha, y, x, lambda)
}
\arguments{
\item{beta}{a p-vector of coefficients}

\item{alpha}{overdispersion parameter of negative binomial}

\item{y}{a n-vector of response}

\item{x}{a n x p matrix of design}

\item{lambda}{mean of Poisson censoring process}
}
\description{
Log-likelihood for negative binomial with interfered responses
}
