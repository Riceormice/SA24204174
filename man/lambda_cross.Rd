% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunc.R
\name{lambda_cross}
\alias{lambda_cross}
\title{lambda_cross}
\usage{
lambda_cross(prior, X, y, betas, lambda, ntimes, k)
}
\arguments{
\item{prior}{A numeric vector representing the prior information for weighting.}

\item{X}{A numeric matrix representing the feature data.}

\item{y}{A numeric vector representing the response variable.}

\item{betas}{A matrix where each column is a candidate coefficient vector.}

\item{lambda}{A numeric vector of candidate `lambda` values to evaluate.}

\item{ntimes}{An integer indicating the number of cross-validation repetitions.}

\item{k}{An integer specifying the number of folds in each cross-validation repetition.}
}
\value{
The optimal `lambda` value that maximizes the cross-validation score.
}
\description{
This function performs cross-validation to select the optimal regularization parameter (`lambda`)
for a given dataset and prior. It uses a custom scoring method to evaluate `lambda` values across
multiple folds and repetitions.
}
\examples{
# Generate synthetic data
set.seed(123)
prior <- c(1, 0.5, 0.2)
X <- matrix(rnorm(300), nrow = 100, ncol = 3)
betas <- t(matrix(rnorm(9), nrow = 3, ncol = 3))
y <- X\%*\%betas[,1]+rnorm(100)
lambda <- seq(0.01, 1, by = 0.1)
ntimes <- 5
k <- 5

# Run the function
best_lambda <- lambda_cross(prior, X, y, betas, lambda, ntimes, k)
print(best_lambda)

}
