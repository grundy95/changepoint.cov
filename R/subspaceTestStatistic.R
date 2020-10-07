#' Test statistic for Subspace method
#'
#' Calculates the Subspace test statistic for all potential changepoint locations within the time series.
#'
#' @param X Data matrix of dimension n by p
#' @param q Dimension of the latent subspace
#' @param msl Minimum segment length between changepoints. NOTE this should be greater than or equal to p
#' 
#' @return A numeric vector containing the test statistic at each potential changepoint location
subspaceTestStat <- function(X,q,msl=dim(X)[2]){
	n <- nrow(X)
	p <- ncol(X)
	X <- purrr::map(1:n,~X[.,])
	calculateSubspaceCost <- subspaceCostCalculator(X,q)
	testStat <- purrr::map_dbl(msl:(n-msl),calculateSubspaceCost)
	return(c(rep(NA,(msl-1)),testStat,rep(NA,msl)))
}


#' Subspace cost calculator
#' 
#' Creates a function that takes potential changepoint location and returns the test statistic as defined in Grundy(2020)
#'
#' @param X Data matrix of dimension n by p
#' @param q Dimension of the latent subspace

#' @return A function used to calculate the test statistic
subspaceCostCalculator <- function(X,q){
	A <- purrr::map(X,~.%*%t(.))
	A <- purrr::accumulate(A,`+`)
	n <- length(X)
	p <- length(X[[1]])
	nullCost <- n*sum(eigen((1/n)*A[[n]],only.values=TRUE,symmetric=TRUE)$values[(q+1):p])
	function(tau){
		eigs1 <- eigen((1/tau)*A[[tau]],only.values=TRUE,symmetric=TRUE)$values
		eigs2 <- eigen((1/(n-tau))*(A[[n]]-A[[tau]]),only.values=TRUE,symmetric=TRUE)$values
		altCost <- tau*sum(eigs1[(q+1):p])+(n-tau)*sum(eigs2[(q+1):p]) 
		return(nullCost-altCost)
	}
}

