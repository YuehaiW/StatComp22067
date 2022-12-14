% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{f_gra1}
\alias{f_gra1}
\title{Compute first-order gradient}
\usage{
f_gra1(x, lambda, s)
}
\arguments{
\item{x}{a positive definite symmetric matrix}

\item{lambda}{a regularization parameter}

\item{s}{the sample covariance}
}
\value{
d the first-order gradient
}
\description{
Compute first-order gradient of the objective function
}
\examples{
\dontrun{
f_gra1(diag(3),1,diag(3))
}
}
