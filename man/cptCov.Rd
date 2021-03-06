% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cptCov.R
\name{cptCov}
\alias{cptCov}
\title{Covariance changepoint detection}
\usage{
cptCov(
  X,
  method = c("Ratio", "CUSUM"),
  threshold = "Asymptotic",
  numCpts = "AMOC",
  msl = 2 * ncol(X),
  thresholdValue = 0.05,
  LRCov = "Bartlett",
  Class = TRUE
)
}
\arguments{
\item{X}{Data matrix of dimension n by p.}

\item{method}{Covariance changepoint method to be used. Choice of "Ratio" or "CUSUM".}

\item{threshold}{Threshold choice for determining significance of changepoints. Choices include:
\itemize{
\item "Asymptotic" - Uses the asymptotic threshold derived for each method. For Ratio method the threshold is log(n). For CUSUM method the threshold is the specified quantile of the standard Normal distribution. The quantile is set by the argument thresholdValue.
\item "Manual"- A user chosen threshold which is contained in the thresholdValue argument. NOTE the normalized test statistics will be compared to the set thresholds - see details for more information.
}
If numCpts is numeric then the threshold is not used as the number of changepoints is known.}

\item{numCpts}{Number of changepoints in the data. Choices include:
\itemize{
\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
\item Numeric - User specified number of changepoints.
}}

\item{msl}{Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.}

\item{thresholdValue}{Either the manual threshold value when threshold="Manual" or the (1-thresholdValue)-quantile of asymptotic distribution of the test statistic when threshold="Asymptotic".}

\item{LRCov}{The long-run covariance estimator to be used for CUSUM method. Currently, only "Bartlett" and "Empirical" are supported. Alternatively, a matrix containing the long-run covariance estimate can be inputted.}

\item{Class}{Logical. If TRUE then an S4 class is returned. If FALSE the estimated changepoints are returned.}
}
\value{
An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
}
\description{
Finds covariance changepoints within multivariate time series data using
either the Ratio method \insertCite{Ryan2020}{changepoint.cov} or a
CUSUM based method \insertCite{Aue2009}{changepoint.cov}.
}
\details{
This is a simple wrapper function for the functions
\code{\link{cptRatio}} and \code{\link{cptCUSUM}}. If no method is specified
then the method used will depend on the dimension of the time series. For
p<10, the CUSUM method will be used and for p>=10 the Ratio method will be
used.
}
\examples{
set.seed(1)
dataAMOC <- wishartDataGeneration(n=100,p=3,tau=50)$data
dataMultipleCpts <- wishartDataGeneration(n=200,p=3,tau=c(50,100,150))$data

ansRatio <- cptCov(X=dataAMOC,method="Ratio")
summary(ansRatio)
plot(ansRatio)

ansCUSUM <- cptCov(X=dataAMOC,method='CUSUM')
show(ansCUSUM)

ansRatio2 <- cptCov(X=dataMultipleCpts,method='Ratio',threshold='Manual',numCpts='BinSeg',
		msl=10,thresholdValue=20)
summary(ansRatio2)
cptsSig(ansRatio2)
plot(ansRatio2)

ansCUSUM2 <- cptCov(X=dataAMOC,method='CUSUM',numCpts=3,
		msl=15,LRCov='Empirical')
summary(ansCUSUM2)
cptsSig(ansCUSUM2)

}
\references{
\insertRef{Ryan2020}{changepoint.cov}

\insertRef{Aue2009}{changepoint.cov}
}
\seealso{
\code{\link{cptRatio}}, \code{\link{cptCUSUM}}, \code{\linkS4class{cptCovariance}}, \code{\link{wishartDataGeneration}}
}
