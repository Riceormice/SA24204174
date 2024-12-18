% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunc.R
\name{Trans.lasso}
\alias{Trans.lasso}
\title{Trans-LASSO}
\usage{
Trans.lasso(X, y, n.vec, I.til, l1 = T)
}
\arguments{
\item{X}{the design matrix}

\item{y}{the response vector}

\item{n.vec}{a vector defining the sample splitting indices; `n.vec[k]` defines the size of the `k`th split}

\item{I.til}{a vector of indices for the transfer dataset to be used for aggregation}

\item{l1}{a boolean flag indicating whether to use lasso regression in the selection steps (default is `TRUE`)}
}
\value{
a list containing:
\item{beta.hat}{the aggregated coefficients using Q-aggregation}
\item{theta.hat}{the weights assigned to models during Q-aggregation}
\item{rank.pi}{the ranks of transfer datasets based on the computed Rhat values}
\item{beta.pool}{the aggregated coefficients using an alternative pool method}
\item{theta.pool}{the weights assigned to models during the alternative pool aggregation}
}
\description{
Perform Transfer Learning with Lasso regression using aggregation and selection
}
\examples{
\dontrun{
X <- matrix(rnorm(1000), 100, 10)
y <- rnorm(100)
n.vec <- c(50, 30, 20)
I.til <- sample(1:50, 10)
result <- Trans.lasso(X, y, n.vec, I.til)
print(result)
}
}
