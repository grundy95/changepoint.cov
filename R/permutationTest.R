#' Permutation Test for Subspace Method
#'
#' Calculates a threshold for the Subspace Method via a permutation test.
#'
#' @param X Data matrix of dimension n by p
#' @param subspaceDim Dimension of the latent subspace
#' @param msl Minimum segment length between changepoints. Note this should be greater than or equal to p
#' @param alpha Significance level of test
#' @param nperm Number of random permutations to be performed
#'
#' @return Numeric containing the threshold value
#'
#' @examples
#' set.seed(1)
#' data <- subspaceDataGeneration(n=100,p=20,subspaceDim=5,tau=50)$data
#' ans <- permutationTest(X=data,subspaceDim=5,msl=20,nperm=100)
#' ans
#'
#' @importFrom stats quantile
#' @export
permutationTest <- function(X,subspaceDim,msl=dim(X)[2],alpha=0.05,nperm=200){
	T <- purrr::rerun(nperm,max(subspaceTestStat(X[sample(nrow(X),replace=FALSE),],subspaceDim,msl),na.rm=TRUE))
	thresh <- quantile(unlist(T),1-alpha)
	return(thresh)
}
