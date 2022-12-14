% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{AR2}
\alias{AR2}
\title{Generate a synthetic data}
\usage{
AR2(n)
}
\arguments{
\item{n}{number of variables}
}
\value{
A the synthetic data
}
\description{
Generate a synthetic data. Specificly, A[i,i]=1, A[i,i-1]=A[i-1,i]=0.5, A[i,i-2]=A[i-2,i]=0.25 and zero otherwise.
}
\examples{
\dontrun{
AR2(30)
}
}
