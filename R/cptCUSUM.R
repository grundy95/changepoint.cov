#' CUSUM based covariance changepoint detection
#'
#' Implements the Aue method for detecting covariance changes in multivariate time series. This method is aimed at low-dimensional time series that can have temporal dependence.
#'
#' @inheritParams cptCov
#' @param errorCheck Logical. If TRUE error checking is performed
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
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
	numCpts <- toupper(numCpts)
	LRCov <- toupper(LRCov)


	n <- nrow(X)
	p <- ncol(X)
	delta <- p*(p+1)/2
	if(numCpts=='AMOC'){
		testStat <- cusumTestStat(X,LRCov,msl)
		T <- (mean(testStat,na.rm=TRUE)-delta/6)/sqrt(delta/45)
		if(threshold=='ASYMPTOTIC'){
			thresh <- qnorm(1-thresholdValue)
		}else if(threshold=='MANUAL'){
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
		return(classInput(X=X,cpts=cpts,method='CUSUM',numCpts=numCpts,testStat=T,threshold=threshold,thresholdValue=thresh,msl=msl,LRCov=LRCov))
	}
}

	


