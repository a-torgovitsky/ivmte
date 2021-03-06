% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testfunctions_covariates.R
\name{sOls3}
\alias{sOls3}
\title{IV-like weighting function, OLS specification 3}
\usage{
sOls3(x, d, j, exx)
}
\arguments{
\item{x}{vector, the value of the covariates other than the
intercept and the treatment indicator.}

\item{d}{0 or 1, indicating treatment or control.}

\item{j}{scalar, position of the component one is interested in
constructing the IV-like weight for.}

\item{exx}{the matrix E[XX']}
}
\value{
scalar.
}
\description{
IV-like weighting function for OLS specification 3.
}
