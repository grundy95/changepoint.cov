#' CUSUM based covariance changepoint detection
#'
#' Implements the Aue method for detecting covariance changes in multivariate time series. This method is aimed at low-dimensional time series that can have temporal dependence.
#'
#' @param X Data matrix of dimension n by p
#' @param threshold Character threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "Normal" - Uses the normality of the asymptotic distribution of the test statistic. The significance level is contained in the thresholdValue parameter.
#'	\item "Manual"- A user chosen threshold which is contained in the thresholdValue parameter. NOTE the normalized test statistic will be compared to the set threshold.
#'}
#' @param noCpts Number of changepoints in the data. Choices include: 
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item "Set" - Number of changepoints in the data is known and is contained in the m parameter. The method returns the m most significant changepoints via a Binary Segmentation framework.
#' }
#' @param LRCov Character stating the long-run covariance estimator to be used. Currently, only "Bartlett" and "Empirical" are supported.
#' @param statType Character stating which test statistic should be used. Choice of "Mean" and "Max".
#' @param thresholdValue Either the significance level of the Normal distribution when using threshold="Normal" or the user defined threshold when using threshold="Manual".
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param Class Logical. If TRUE then an S4 class is returned. Else just the estimated changepoints are returned.

#' @return An object of S4 class \code{\link{cptCov}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export

cptCUSUM <- function(X,threshold='Normal',noCpts='AMOC',LRCov='Bartlett',statType='Mean',thresholdValue=0.05,msl=dim(X)[2],Class=TRUE){
	cusumErrorChecks(X=X,threshold=threshold,noCpts=noCpts,LRCov=LRCov,statType=statType,thresholdValue=thresholdValue,msl=msl,Class=Class)
	n <- dim(X)[1]
	p <- dim(X)[2]
	delta <- p*(p+1)/2
	if(noCpts=='AMOC'){
		testStat <- cusumTestStat(X,LRCov,msl)
		if(statType=='Mean'){
			T <- (mean(testStat,na.rm=TRUE)-delta/6)/sqrt(delta/45)
		}else if(statType=='Max'){
			T <- (max(testStat,na.rm=TRUE)-delta/4)/sqrt(delta/8)
		}

		if(threshold=='Normal'){
			thresh <- qnorm(1-thresholdValue)
		}else if(threshold=='Manual'){
			thresh <- thresholdValue
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
		return(classInput(X=X,cpts=cpts,method='CUSUM',noCpts=noCpts,testStat=T,threshold=threshold,thresholdValue=thresh,msl=msl,LRCov=LRCov,statType=statType))
	}
}

	


