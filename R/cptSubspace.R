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
#' @param numCpts Number of changepoints in the data. Choices include: 
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item "Set" - Number of changepoints in the data is known and is contained in the m parameter. The method returns the m most significant changepoints via a Binary Segmentation framework.
#' }
#' @param thresholdValue Either the significance level of the permutation test when using threshold="PermTest" or the user defined threshold when using threshold="Manual".
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param nperm Only required for threshold="PermTest". Number of permutations to use in the permutation test.
#' @param Class Logical. If TRUE then an S4 class is returned. Else just the estimated changepoints are returned.

#' @return An object of S4 class "cptCovSubspace" is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export
cptSubspace <- function(X,q,threshold='PermTest',numCpts='AMOC',thresholdValue=0.05,msl=ncol(X),nperm=200,Class=TRUE){
	if(is.data.frame(X)){
		X <- as.matrix(X)
	}
	subspaceErrorChecks(X=X,q=q,threshold=threshold,numCpts=numCpts,thresholdValue=thresholdValue,msl=msl,nperm=nperm,Class=Class)
	threshold <- toupper(threshold)
	numCpts <- toupper(numCpts)
	n <- nrow(X)
	p <- ncol(X)
	if(numCpts=='AMOC'){
		testStat <- subspaceTestStat(X,q,msl)
		T <- max(testStat,na.rm=TRUE)
		if(threshold=='PERMTEST'){
			thresh <- permutationTest(X,q,msl,thresholdValue,nperm)
		}else if(threshold=='MANUAL'){
			thresh <- thresholdValue
			nperm <- 0
		}
		if(T>thresh){
			cpts <- c(which.max(testStat),n)
		}else{
			cpts <- c(n)
		}
	}
	if(Class==FALSE){
		return(list('cpts'=cpts))
	}else{
		return(classInput(X=X,cpts=cpts,method='Subspace',numCpts=numCpts,testStat=T,threshold=threshold,thresholdValue=thresh,msl=msl,q=q,nperm=nperm))
	}
}

