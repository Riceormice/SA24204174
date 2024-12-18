% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{un}
\alias{un}
\title{Pairwise U-statistic computation}
\usage{
un(x, y)
}
\arguments{
\item{x}{A numeric matrix where each row represents an observation and each column represents a feature.}

\item{y}{A numeric vector representing the response variable, with length equal to the number of rows in `x`.}
}
\value{
A numeric vector of length equal to the number of columns in `x`, representing the computed U-statistic.
}
\description{
Computes a U-statistic based on pairwise comparisons of rows in a predictor matrix `x` and response vector `y`.
}
\examples{
\dontrun{
# Example data
x <- matrix(rnorm(1000), ncol = 5) # 200 rows, 5 columns
y <- rnorm(200)
# Compute U-statistic
u_stat <- un(x, y)
print(u_stat)
}
}
