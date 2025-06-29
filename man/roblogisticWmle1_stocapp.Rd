% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roblogistic.R
\name{roblogisticWmle1_stocapp}
\alias{roblogisticWmle1_stocapp}
\title{Stochastic approximation for robust logistic regression with inconsistent initial estimator with Tukey's weights}
\usage{
roblogisticWmle1_stocapp(
  x,
  thetastart,
  c = 4.685061,
  maxit = 10000,
  tol = 1e-07,
  verbose = FALSE,
  seed = 321
)
}
\arguments{
\item{x}{a n x p matrix of design}

\item{thetastart}{an inconsistent estimator (also used as starting values)}

\item{c}{tuning parameter for Tukey's weight (default value is 4.685061)}

\item{maxit}{max number of iteration for IRWLS}

\item{tol}{tolerance for stopping criterion}

\item{verbose}{print info}

\item{seed}{for random number generator}

\item{H}{number of estimators for Monte Carlo approximation}
}
\description{
Stochastic approximation for robust logistic regression with inconsistent initial estimator with Tukey's weights
}
