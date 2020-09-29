#' Fisher matrix based covariance changepoint detection
#'
#' Implements the Ryan method for detecting covariance changes in multivariate time series. This method is aimed at independent high-dimensional time series.
#'
#' @param X Data matrix of dimension n by p
#' @param threshold Character threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "Asymptotic" - Uses the asymptotic distribution of the test statistic. See details for more information
#'	\item "Manual"- A user chosen threshold which is contained in the thresholdValue parameter. 
#'}
#' @param noCpts Number of changepoints in the data. Choices include: 
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item "Set" - Number of changepoints in the data is known and is contained in the m parameter. The method returns the m most significant changepoints via a Binary Segmentation framework.
#' }
#' @param thresholdValue User defined threshold when using threshold="Manual".
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param Class Logical. If TRUE then an S4 class is returned. Else just the estimated changepoints are returned.

#' @return An object of S4 class \code{\link{cptCov}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export

cptFisher <- function(X,threshold='Asymptotic',noCpts='AMOC',thresholdValue=NA,msl=dim(X)[2],Class=TRUE){
	fisherErrorChecks(X=X,threshold=threshold,noCpts=coCpts,thresholdValue=thresholdValue,msl=msl,Class=Class)
	n <- nrow(X)
	p <- ncol(X)
	if(noCpts=='AMOC'){
		testStat <- fisherTestStat(X,msl)
		T <- max(testStat,na.rm=TRUE)
		if(threshold=="Asymptotic"){
			thresh <- log(n)
		}else if(threshold=="Manual"){
			thresh <- thesholdValue
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
		return(classInput(X=X,cpts=cpts,method='Fisher',noCpts=noCpts,testStat=T,threshold=threshold,thresholdValue=thresh,msl=msl))
	}
}

