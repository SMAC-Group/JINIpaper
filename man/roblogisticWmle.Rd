% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{roblogisticWmle}
\alias{roblogisticWmle}
\title{Robust logistic regression estimator with Tukey's weights}
\usage{
roblogisticWmle(
  y,
  x,
  start,
  c = 4.685061,
  maxit = 200L,
  tol = 1e-07,
  verbose = FALSE
)
}
\arguments{
\item{y}{a vector of responses}

\item{x}{a n x p matrix of design}

\item{c}{tuning parameter for Tukey's weight (default value is 4.685061)}

\item{maxit}{max number of iteration for IRWLS}

\item{tol}{tolerance for stopping criterion}

\item{verbose}{print info}

\item{beta}{a p-vector of parameter (starting values)}
}
\description{
Robust logistic regression estimator with Tukey's weights
}
