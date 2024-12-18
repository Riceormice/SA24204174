% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfunc.R
\name{gic}
\alias{gic}
\title{Calculate GIC}
\usage{
gic(beta, n, an, u)
}
\arguments{
\item{beta}{the estimator}

\item{n}{the number of samples}

\item{an}{the function using in gic}

\item{u}{the U-statistics}
}
\value{
gic
}
\description{
Calculate GIC of estimator
}
\examples{
\dontrun{
gic(c(1,2),5,log,c(1,2))

}
}
