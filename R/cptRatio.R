#' Ratio method for covariance changepoint detection
#'
#' Implements the Ryan method for detecting covariance changes in multivariate time series. This method is aimed at independent high-dimensional time series.
#'
#' @inheritParams cptCov
#' @param errorCheck Logical. If TRUE error checking is performed
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#'
#' @examples
#' set.seed(1)
#' dataAMOC <- wishartDataGeneration(n=100,p=5,tau=50)$data
#'
#' ansRatio <- cptRatio(X=dataAMOC)
#' summary(ansRatio)
#' 
#' ansRatio2 <- cptCov(X=dataAMOC,threshold='Manual',numCpts='AMOC',msl=10,thresholdValue=20)
#' summary(ansRatio2)
#'
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
	}
	threshold <- toupper(threshold)
	if(is.character(numCpts)){
		numCpts <- toupper(numCpts)
	}


	n <- nrow(X)
	p <- ncol(X)
	if(threshold=="ASYMPTOTIC"){
		thresholdValue <- log(n)
	}
	if(numCpts=='AMOC'){
		testStat <- ratioTestStat(X,msl)
		T <- max(testStat,na.rm=TRUE)
		if(T>thresholdValue){
			cpts <- c(which.max(testStat),n)
			cptsSig <- data.frame('cpts'=cpts[1],'T'=T,'thresholdValue'=thresholdValue)
		}else{
			cpts <- c(n)
			cptsSig <- data.frame('cpts'=which.max(testStat),'T'=T,'thresholdValue'=thresholdValue)
		}
	}else{
		if(is.numeric(numCpts)){
			thresholdValue <- 0
			cptsSig <- binSeg(X,method='RATIO',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=numCpts)
		}else{
			cptsSig <- binSeg(X,method='RATIO',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1)
		}
		cpts <- c(cptsSig$cpts,n)
		cpts <- sort(cpts[!is.na(cpts)])
	}
	if(Class==FALSE){
		return(cpts)
	}else{
		return(classInput(X=X,cpts=cpts,method='Ratio',numCpts=numCpts,testStat=cptsSig$T,threshold=threshold,thresholdValue=thresholdValue,msl=msl))
	}
}

