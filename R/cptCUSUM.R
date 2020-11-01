#' CUSUM based covariance changepoint detection
#'
#' Implements the \insertCite{Aue2009;textual}{changepoint.cov} method for detecting covariance changes in multivariate time series. This method is aimed at low-dimensional time series that can have temporal dependence.
#'
#' This function calculates the test statistic, T, described in \insertCite{Aue2009;textual}{changepoint.cov}, specifically the sum of the test statistics for all potential changepoint locations; changepoint localization uses the maximum of these test statistics. T is then normalized so that its asymptotic distribution is a standard Normal. This normalised test statistic is then compared to the defined threshold (either the specified quantile of the Normal distribution if the threshold is set as asymptotic or the manual threshold). If multiple changepoints are possible then the Binary Segmentation algorithm is used to detect multiple changes. The long run covariance estimation and its inversion can be unstable if the dimension of the time series is large. In this scenario we recommend using the \code{\link{cptRatio}} function.
#'
#' @inheritParams cptCov
#' @param errorCheck Logical. If TRUE error checking is performed.
#'
#' @return An object of S4 class \code{\linkS4class{cptCovariance}} is returned. If Class="FALSE", the vector of changepoints are returned.
#'
#'@references
#' \insertRef{Aue2009}{changepoint.cov}
#'
#' @seealso \code{\link{cptCov}}, \code{\linkS4class{cptCovariance}}, \code{\link{wishartDataGeneration}}, \code{\link{cusumTestStat}}
#'
#' @examples
#' set.seed(1)
#' dataAMOC <- wishartDataGeneration(n=100,p=2,tau=50)$data
#'
#' ansCUSUM <- cptCov(X=dataAMOC)
#' show(ansCUSUM)
#'
#' ansCUSUM2 <- cptCov(X=dataAMOC,threshold='Manual',numCpts='AMOC',
#'			msl=15,thresholdValue=15,LRCov='Empirical')
#' summary(ansCUSUM2)
#'
#' @importFrom stats qnorm 
#'
#' @include cptCovarianceClass.R
#'
#' @export

cptCUSUM <- function(X,threshold='Asymptotic',numCpts='AMOC',msl=2*ncol(X),LRCov='Bartlett',thresholdValue=0.05,errorCheck=TRUE,Class=TRUE){
	if(!is.logical(errorCheck)){
		stop("errorCheck should be logical")
	}
	if(errorCheck){
		if(is.data.frame(X)){
			X <- as.matrix(X)
		}
		cusumErrorChecks(X=X,threshold=threshold,numCpts=numCpts,LRCov=LRCov,thresholdValue=thresholdValue,msl=msl,Class=Class)
	}
	threshold <- toupper(threshold)
	if(is.character(numCpts)){
		numCpts <- toupper(numCpts)
	}
	if(is.character(LRCov)){
		LRCov <- toupper(LRCov)
	}

	n <- nrow(X)
	p <- ncol(X)
	if(threshold=='ASYMPTOTIC'){
		thresholdValue <- qnorm(1-thresholdValue)
	}
	if(numCpts=='AMOC'){
		testStat <- cusumTestStat(X,LRCov,msl)
		T <- mean(testStat,na.rm=TRUE)
		if(T>thresholdValue){
			cpts <- c(which.max(testStat),n)
			cptsSig <- data.frame('cpts'=cpts[1],'T'=T,'thresholdValue'=thresholdValue,'significant'=TRUE)
		}else{
			cpts <- c(n)
			cptsSig <- data.frame('cpts'=which.max(testStat),'T'=T,'thresholdValue'=thresholdValue,'significant'=FALSE)
		}
	}else{
		if(is.numeric(numCpts)){
			thresholdValue <- 0
			cptsSig <- binSeg(X,method='CUSUM',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=numCpts,LRCov=LRCov)
		}else{
			cptsSig <- binSeg(X,method='CUSUM',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1,LRCov=LRCov)
		}
		cpts <- c(cptsSig$cpts[cptsSig$significant],n)
	}
	if(Class==FALSE){
		return(cpts)
	}else{
		return(classInput(X=X,cpts=cpts,method='CUSUM',numCpts=numCpts,cptsSig=cptsSig,threshold=threshold,thresholdValue=unique(cptsSig$thresholdValue),msl=msl,LRCov=LRCov))
	}
}

	


