#' Test statistic for Subspace method
#'
#' Calculates the Subspace test statistic for all potential changepoint locations within the time series. See \code{\link{cptSubspace}} for more details.
#'
#' See \code{\link{cptSubspace}}.
#'
#' @param X Data matrix of dimension n by p.
#' @param subspaceDim Dimension of the latent subspace.
#' @param msl Minimum segment length between changepoints. NOTE this should be greater than or equal to p.
#' 
#' @return A numeric vector containing the test statistic at each potential changepoint location
#'
#' @references
#' \insertRef{Grundy2020}{changepoint.cov}
#'
#' @seealso \code{\link{cptSubspace}}
#'
#' @examples
#' set.seed(1)
#' data <- subspaceDataGeneration(n=100,p=20,subspaceDim=5,tau=50,changeSize=0.5*sqrt(5))$data
#' ans <- subspaceTestStat(X=data,subspaceDim=5,msl=30)
#' which.max(ans)
#'
#' @export
subspaceTestStat <- function(X,subspaceDim,msl){
	n <- nrow(X)
	p <- ncol(X)
	X <- purrr::map(1:n,~X[.,])
	calculateSubspaceCost <- subspaceCostCalculator(X,subspaceDim)
	testStat <- purrr::map_dbl(msl:(n-msl),calculateSubspaceCost)
	return(c(rep(NA,(msl-1)),testStat,rep(NA,msl)))
}


#' Subspace cost calculator
#' 
#' Creates a function that takes potential changepoint location and returns the test statistic as defined in \insertCite{Grundy2020;textual}{changepoint.cov}.
#'
#' @param X Data matrix of dimension n by p
#' @param subspaceDim Dimension of the latent subspace
#'
#' @return A function used to calculate the test statistic
#'
#' @seealso \code{\link{subspaceTestStat}}
#'
#' @keywords internal
subspaceCostCalculator <- function(X,subspaceDim){
	A <- purrr::map(X,~.%*%t(.))
	A <- purrr::accumulate(A,`+`)
	n <- length(X)
	p <- length(X[[1]])
	nullCost <- n*sum(eigen((1/n)*A[[n]],only.values=TRUE,symmetric=TRUE)$values[(subspaceDim+1):p])
	function(tau){
		eigs1 <- eigen((1/tau)*A[[tau]],only.values=TRUE,symmetric=TRUE)$values
		eigs2 <- eigen((1/(n-tau))*(A[[n]]-A[[tau]]),only.values=TRUE,symmetric=TRUE)$values
		altCost <- tau*sum(eigs1[(subspaceDim+1):p])+(n-tau)*sum(eigs2[(subspaceDim+1):p]) 
		return(nullCost-altCost)
	}
}

