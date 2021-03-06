% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cusumTestStat.R
\name{cusumVech}
\alias{cusumVech}
\title{Column stacking of covariance matrix}
\usage{
cusumVech(X)
}
\arguments{
\item{X}{List of data where each slot is a time point}
}
\value{
An n by p(p+1)/2 matrix
}
\description{
Calculates the covariance matrix for each time point and stacks the columns using \code{\link{vech}}. Used in CUSUM method.
}
\seealso{
\code{\link{cusumTestStat}}
}
\keyword{internal}
