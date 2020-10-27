#' Test statistic for CUSUM method
#'
#' Calculates the test statistic for all potential changepoint locations within the time series.
#'
#' @param X Data matrix of dimension n by p
#' @param LRCov A character describing the Long-run covariance estimator to be used. Currently, only "Bartlett" and "Empirical" are available, see \code{\link{cptCUSUM}} for more details.
#' @param msl A numeric giving the minimum segment length between changepoints. NOTE this should be greater than or equal to p.
#'
#' @return A numeric vector containing the test statistic at each potential changepoint location.
#'
#' @examples
#' data <- wishartDataGeneration(n=100,p=3,tau=50)$data
#' ans <- cusumTestStat(X=data,LRCov='Empirical',msl=20)
#' which.max(ans)
#'
#' @importFrom stats cov
#' @export
cusumTestStat <- function(X,LRCov,msl){
	n <- nrow(X)
	p <- ncol(X)
	delta <- (p*(p+1))/2
	if((is.matrix(LRCov))&&(!((ncol(LRCov)==delta)&&(nrow(LRCov)==delta)))){
		stop("Dimension of manual LRCov is not compatible with data")
	}	
	if(is.character(LRCov)){
		LRCov <- toupper(LRCov)
	}
	X <- purrr::map(1:n,~X[.,])
	Xvech <- cusumVech(X)
	if((is.character(LRCov))&&(LRCov=='BARTLETT')){
		#Should we be multiplying by n or n-delta?
		cusumCov <- tryCatch({
			sandwich::lrvar(Xvech,kernel='Bartlett')*n
		},
		warning=function(cond){
			stop("Long run covariance estimation not possible, try a different long run covariance estimator or the Ratio method")
		},silent=TRUE)
	}else if((is.character(LRCov))&&(LRCov=='EMPIRICAL')){
		cusumCov <- cov(Xvech)
	}else{
		cusumCov <- LRCov
	}
	calculateCusum <- CusumCalculator(X)
	Cusum <- purrr::map(msl:(n-msl),calculateCusum)
	cusumStat <- tryCatch({
		purrr::map_dbl(Cusum,~.%*%solve(cusumCov,.))},
		error=function(cond){
			stop("Long run covariance estimator is not invertible. This is most likely due to the dimension of the data being too large. Please try the Ratio method")
		}
		)
	return(c(rep(NA,(msl-1)),(cusumStat-delta/6)/(sqrt(delta/45)),rep(NA,msl)))
}

#' CUSUM Calculator
#'
#' Creates a function that takes a changepoint location as an input and returns the CUSUM statistic assuming a changepoint at the given location.
#'
#' @param X List of data where each slot is a time point
#'
#' @return A function used to calculate the CUSUM statistic
CusumCalculator <- function(X){
	A <- purrr::map(X,~.%*%t(.))
	A <- purrr::accumulate(A,`+`)
	A <- purrr::map(A,vech)
	n <- length(X)
	function(tau){
		return((1/sqrt(n))*(A[[tau]]-(tau/n)*A[[n]]))
	}
}

#' Column stacking
#'
#' Stacks the columns below the diagonal of a symmetric matrix
#'
#' @param X A data matrix of dimension n by p
#'
#' @return A numeric vector of length p(p+1)/2
vech <- function(X){
	return(as.vector(X[upper.tri(X,diag=TRUE)]))
}

#' Column stacking of covariance matrix
#'
#' Calculates the covariance matrix for each time point and stacks the columns using \code{\link{vech}}
#'
#' @param X List of data where each slot is a time point
#'
#' @return An n by p(p+1)/2 matrix 
cusumVech <- function(X){
	A <- purrr::map(X,~.%*%t(.))
	A <- purrr::map(A,vech)
	A <- purrr::reduce(A,rbind)
	return(A)
}
