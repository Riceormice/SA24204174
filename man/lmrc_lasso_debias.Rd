% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lmrc_lasso_debias}
\alias{lmrc_lasso_debias}
\title{LMRC Lasso Debiased Estimator}
\usage{
lmrc_lasso_debias(
  n,
  p,
  sigma,
  L,
  s,
  u,
  beta_ini,
  beta_ori,
  w,
  denominator,
  lambda = 0.1,
  alpha = 0.15,
  max_iter = 100L,
  tol = 5e-06,
  l1_lambda = 0.01
)
}
\arguments{
\item{n}{Integer. The number of samples.}

\item{p}{Integer. The number of predictors (features).}

\item{sigma}{Eigen::MatrixXd. The covariance matrix of the predictors.}

\item{L}{Double. The Lipschitz constant for the gradient.}

\item{s}{Double. Indicator for positive (1) or negative (-1) direction. If 0, computes both and chooses the better one.}

\item{u}{Eigen::VectorXd. The gradient vector.}

\item{beta_ini}{Eigen::VectorXd. Initial beta estimate.}

\item{beta_ori}{Eigen::VectorXd. Original beta estimate.}

\item{w}{Eigen::VectorXd. A weighting vector.}

\item{denominator}{Eigen::VectorXd. A denominator vector for scaling.}

\item{lambda}{Double. Regularization parameter for debiasing. Default is 0.1.}

\item{alpha}{Double. Step size parameter. Default is 0.15.}

\item{max_iter}{Integer. Maximum number of iterations. Default is 100.}

\item{tol}{Double. Convergence tolerance. Default is 5e-6.}

\item{l1_lambda}{Double. Regularization parameter for L1 penalty. Default is 0.01.}
}
\value{
Eigen::VectorXd. The debiased beta estimate.
}
\description{
Computes a debiased estimator using LMRC Lasso for a given set of parameters.
}
\examples{
\dontrun{
n <- 100
p <- 50
beta_eg <- c(rep(0.5,10),p-10)
Sig.X <- diag(rep(1,p))
X <- rmvnorm(n.vec[k], rep(0, p), Sig.X)
y <- X\%*\%beta_eg + rnorm (n, 0, 1))
L <- 1.0
u <- un(X,y)
ini <- rep(1, p)
ori <- rep(0, p)
w <- rep(1, p)
weight <- rep(1, p)
lam1 <- 0.1
alpha <- 0.15
iter <- 100
tol <- 5e-6
lam2 <- 0.01
s = sign(u[1])
bt <- lmrc_lasso_debias(n, p, cov(X), L, s, u, ini, ori, w, weight, lam1, alpha, iter, tol, lam2)
cat(bt)
}
}
