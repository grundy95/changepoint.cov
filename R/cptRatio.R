#' Ratio method for covariance changepoint detection
#'
#' Implements the Ryan method for detecting covariance changes in multivariate time series. This method is aimed at independent high-dimensional time series.
#'
#' @inheritParams cptCov
#' @param errorCheck Logical. If TRUE error checking is performed
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export

cptRatio <- function(X,threshold='Asymptotic',numCpts='AMOC',msl=2*ncol(X),thresholdValue=0,errorCheck=TRUE,Class=TRUE){
	if(!is.logical(errorCheck)){
		stop("errorCheck should be logical")
	}
	if(errorCheck){
		if(is.data.frame(X)){
			   X <- as.matrix(X)
		}
		ratioErrorChecks(X=X,threshold=threshold,numCpts=numCpts,thresholdValue=thresholdValue,msl=msl,Class=Class)
		threshold <- toupper(threshold)
		numCpts <- toupper(numCpts)
	}


	n <- nrow(X)
	p <- ncol(X)
	if(numCpts=='AMOC'){
		testStat <- ratioTestStat(X,msl)
		T <- max(testStat,na.rm=TRUE)
		if(threshold=="ASYMPTOTIC"){
			thresh <- log(n)
		}else if(threshold=="MANUAL"){
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
		return(classInput(X=X,cpts=cpts,method='Ratio',numCpts=numCpts,testStat=T,threshold=threshold,thresholdValue=thresh,msl=msl))
	}
}

