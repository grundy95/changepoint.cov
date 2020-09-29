#' Detecting changes in subspace
#'
#' This method allows the detection of subspace changes in multivariate time series data.
#'
#' @param X Data matrix of dimension n by p
#' @param q Dimension of the latent subspace
#' @param threshold Threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "PermTest" - Permutation test is performed using the number of permutations and significance level contained in the nperm and thresholdValue parameters respectively
#'	\item "Manual" - A user chosen threshold is used which is contained in the thresholdValue parameter.
#' }
#' @param noCpts Number of changepoints in the data. Choices include: 
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item "Set" - Number of changepoints in the data is known and is contained in the m parameter. The method returns the m most significant changepoints via a Binary Segmentation framework.
#' }
#' @param thresholdValue Either the significance level of the permutation test when using threshold="PermTest" or the user defined threshold when using threshold="Manual".
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param nperm Only required for threshold="PermTest". Number of permutations to use in the permutation test.
#' @param m Only required for noCpts="Set". Known number of changepoints in the data.
#' @param Class Logical. If TRUE then an S4 class is returned. Else just the estimated changepoints are returned.

#' @return An object of S4 class "cptCovSubspace" is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export
cptSubspace <- function(X,q,threshold='PermTest',noCpts='AMOC',thresholdValue=0.05,msl=dim(X)[2],nperm=200,m=NA,Class=TRUE){
	subspaceErrorChecks(X=X,q=q,threshold=threshold,noCpts=noCpts,thresholdValue=thresholdValue,msl=msl,nperm=nperm,m=m,Class)
	n <- dim(X)[1]
	p <- dim(X)[2]
	if(noCpts=='AMOC'){
		testStat <- subspaceTestStat(X,q,msl)
		if(threshold=='PermTest'){
			thresh <- permutationTest(X,q,msl,thresholdValue,nperm)
		}else if(threshold=='Manual'){
			thresh <- thresholdValue
			nperm <- NA
		}
		if(testStat$T>thresh){
			cpts <- c(testStat$cpt,n)
		}else{
			cpts <- c(n)
		}
	}
	if(Class==FALSE){
		return(list('cpts'=cpts))
	}else{
		return(classInput(X=X,cpts=cpts,method='Subspace',noCpts=noCpts,testStat=testStat$T,threshold=threshold,thresholdValue=thresh,msl=msl,q=q,nperm=nperm))
	}
}

