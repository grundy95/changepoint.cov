#' changepoint.cov: A package for covariance changepoint detection
#'
#' The changepoint.cov package provides 3 methods for covariance changepoint detection.
#'
#' @section Ratio method:
#' The Ratio method is implemented in the \code{\link{cptRatio}} function. This function implements the method described in \insertCite{Ryan2020;textual}{changepoint.cov} and is designed for high-dimensional independent time series.
#' 
#' @section CUSUM method:
#' The CUSUM method is implemented in the \code{\link{cptCUSUM}} function. This function implements the method described in \insertCite{Aue2009;textual}{changepoint.cov} and is designed for low-dimensional time series that could exhibit temporal dependence.
#'
#' @section Subspace method:
#' The subspace method is implemented in the \code{\link{cptSubspace}} function. This function implements the method described in \insertCite{Grundy2020;textual}{changepoint.cov} and is designed for time series which lie in a low-dimensional subspace of known dimension.
#'
#' @seealso \code{\link{cptCov}}, \code{\link{cptRatio}}, \code{\link{cptCUSUM}}, \code{\link{cptSubspace}}, \code{\linkS4class{cptCovariance}}, 
#'
#' @docType package
#' @name changepoint.cov
#' @aliases changepoint.cov-package
#'
#' @importFrom Rdpack reprompt
NULL
