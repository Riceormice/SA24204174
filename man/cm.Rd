% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunc.R
\name{cm}
\alias{cm}
\title{CM}
\usage{
cm(beta, beta0)
}
\arguments{
\item{beta}{the estimated coefficients}

\item{beta0}{the true coefficients}
}
\value{
a binary value: 1 if the support of `beta` matches exactly with the support of `beta0`, and 0 otherwise
}
\description{
Calculate the correctness metric for an estimator
}
\examples{
\dontrun{
beta <- c(0, 1, 0, 1)
beta0 <- c(0, 1, 0, 1)
cm(beta, beta0)  # Returns 1

beta <- c(0, 1, 0, 1)
beta0 <- c(0, 1, 1, 1)
cm(beta, beta0)  # Returns 0
}
}
