#' Test statistic for Subspace method
#'
#' Calculates the test statistic, most likely changepoint location and the null and alternative costs for the Subspace Covariance changepoint method.
#'
#' @param X Data matrix of dimension n by p
#' @param q Dimension of the latent subspace
#' @param msl Minimum segment length between changepoints. NOTE this should be greater than or equal to p
#' 
#' @return A list with slots:
#' \itemize{
#'	\item T - Test statistic 
#'	\item cpt - Most likely changepoint position
#' 	\item nullCost - Cost of data with no changes
#'	\item altCost - List with cost of data with a changepoint at the slot location
#' }
subspaceTestStat <- function(X,q,msl=dim(X)[2]){
	n <- dim(X)[1]
	p <- dim(X)[2]
	nullCost <- subspaceCost(X,q)
	altCost <- rep(NA,n)
	altCost <- purrr::map(msl:(n-msl),~subspaceCost(X[1:.,],q)+subspaceCost(X[(.+1):n,],q))
	tau <- which.min(altCost)
	return(list('T'=nullCost-altCost[[tau]],'cpt'=tau,'nullCost'=nullCost,'altCost'=altCost))
}

#' Cost of data for Subspace method
#' 
#' Calculates the cost of data assuming a latent subspace of dimension q
#'
#' @param X Data matrix of dimension n by p
#' @param q Dimension of the latent subspace

#' @return Cost of data
subspaceCost <- function(X,q){
	n <- dim(X)[1]
	p <- dim(X)[2]
	Sigma <- (1/n)*t(X)%*%X
	cost <- n*sum(eigen(Sigma,only.values=TRUE,symmetric=TRUE)$values[(q+1):p])
	return(cost)
}

