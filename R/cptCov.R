#' Covariance changepoint detection
#'
#' Finds covariance changepoints within multivariate time series data using the Ratio method of Ryan (2020) or the CUSUM method of Aue(2009)
#'
#' @param X Data matrix of dimension n by p
#' @param method Covariance changepoint method to be used. Choice of "Ratio" or "CUSUM".
#' @param threshold Threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "Asymptotic" - Uses the asymptotic threshold derived for each method. For method="Ratio" the threshold is log(n). For method="CUSUM" the threshold is the specified quantile of the standard normal distribution. The quantile is set by the argument thresholdValue.
#'	\item "Manual"- A user chosen threshold which is contained in the thresholdValue argument. NOTE the normalized test statistics will be compared to the set thresholds.
#'}
#' @param numCpts Number of changepoints in the data. Choices include: 
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item "Set" - Number of changepoints in the data is known and is contained in the m parameter. The method returns the m most significant changepoints via a Binary Segmentation framework.
#' }
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param LRCov The long-run covariance estimator to be used for CUSUM method. Currently, only "Bartlett" and "Empirical" are supported.
#' @param thresholdValue Either the manual threshold value when threshold="Manual" or the (1-thresholdValue)-quantile of asymptotic distribution of the CUSUM test statistic when method="CUSUM" and threshold="Asymptotic" 
#' @param Class Logical. If TRUE then an S4 class is returned. Else just the estimated changepoints are returned.
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export
cptCov <- function(X,method=c("Ratio","CUSUM"),threshold="Asymptotic",numCpts='AMOC',msl=2*ncol(X),LRCov='Bartlett',thresholdValue=0.05,Class=TRUE){
	if(!is.character(method)){
		stop("method not recognized: Please choose between 'Ratio' and 'CUSUM'")
	}
	if((length(method)>1)&&(method[1]=="Ratio")&&(method[2]=="CUSUM")){
		if(ncol(X)>20){
			method <- "Ratio"
			warning("no method was chosen. As p>20 the Ratio method will be implemented")
		}else{
			method <- "CUSUM"
			warning("no method was chosen. As p<=20 the CUSUM method will be implemented")
		}
	}	
	if(length(method)!=1){
		stop("only one method can be implemented at once")
	}
	method <- toupper(method)
	if((method=='SUBSPACE')||(method=='GRUNDY')){
		stop('For subspace changepoint detection use the function cptSubspace. See ?cptSubspace for details')
	}
	if(method=="RYAN"){
		method <- "RATIO"
	}
	if(method=='AUE'){
		method <- "CUSUM"
	}
	if((method!="RATIO")&&(method!="CUSUM")){
		stop("method not recognized: Please choose between 'Ratio' and 'CUSUM'")
	}
	if(is.data.frame(X)){
		X <- as.matrix(X)
	}

	if(method=="RATIO"){
		ratioErrorChecks(X=X,threshold=threshold,numCpts=numCpts,thresholdValue=thresholdValue,msl=msl,Class=Class)
		if(LRCov!="Bartlett"){
			warning("Long run covariance estimator is not used in Ratio method")
		}
		threshold <- toupper(threshold)
		numCpts <- toupper(numCpts)
		ans <- cptRatio(X=X,threshold=threshold,numCpts=numCpts,msl=msl,thresholdValue=thresholdValue,errorCheck=FALSE,Class=Class)
	}else{
		cusumErrorChecks(X=X,threshold=threshold,numCpts=numCpts,msl=msl,LRCov=LRCov,thresholdValue=thresholdValue,Class=Class)
		threshold <- toupper(threshold)
		numCpts <- toupper(numCpts)
		LRCov <- toupper(LRCov)
		ans <- cptCUSUM(X=X,threshold=threshold,numCpts=numCpts,msl=msl,LRCov=LRCov,thresholdValue=thresholdValue,errorCheck=FALSE,Class=Class)
	}
	return(ans)
}


