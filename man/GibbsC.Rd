% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GibbsC}
\alias{GibbsC}
\title{A Gibbs sampler using Rcpp}
\usage{
GibbsC(N, mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 1, rho = 0.9)
}
\arguments{
\item{N}{the number of samples}

\item{mu1}{mean of the first random variable}

\item{mu2}{mean of the second random variable}

\item{sigma1}{standard deviation of the first random variable}

\item{sigma2}{standard deviation of the second random variable}

\item{rho}{correlation of the two random variables}
}
\value{
X random sample of size \code{n}
}
\description{
A Gibbs sampler using Rcpp
}
\examples{
\dontrun{
gc <- GibbsC(100)
par(mfrow=c(2,1));
plot(gc[,1],type='l')
plot(gc[,2],type='l')
}
}
