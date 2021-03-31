#' Ratio method for covariance changepoint detection
#'
#' Implements the \insertCite{Ryan2020;textual}{changepoint.cov} method for
#' detecting covariance changes in multivariate time series. This method is
#' aimed at independent high-dimensional time series.
#'
#' This function calculates the test statistic, T, described in
#' \insertCite{Ryan2020;textual}{changepoint.cov}. Using results from Random
#' Matrix Theory the test statistic is normalised by it's asymptotic expectation
#' and variance so that it follows a standard Normal distribution. Following the
#' paper, the threshold \code{log(n)} is used if the threshold is set as
#' asymptotic, else the user defined manual threshold is used. If multiple
#' changepoints are possible then the Binary Segmentation algorithm is used to
#' detect multiple changes. If the minimum segment length is too small then the
#' numerical integration performed in the normalization of the test statistic can
#' be unstable. In this scenario the minimum segment length will be automatically
#' increased. This method is designed for independent time series, if the time
#' series contains temporal dependence we recommend using the
#' \code{\link{cptCUSUM}} function.
#'
#' @inheritParams cptCov
#' @param errorCheck Logical. If TRUE error checking is performed
#'
#' @return An object of S4 class \code{\linkS4class{cptCovariance}} is returned. If Class="FALSE", the vector of changepoints are returned.
#'
#' @references
#' \insertRef{Ryan2020}{changepoint.cov}
#'
#' @seealso \code{\link{cptCov}}, \code{\linkS4class{cptCovariance}}, \code{\link{wishartDataGeneration}}, \code{\link{ratioTestStat}}
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
#' @include cptCovariance-class.R
#'
#' @export

cptRatio <- function(X, threshold='Asymptotic', numCpts='AMOC', msl=2*ncol(X),
                     thresholdValue=0, errorCheck=TRUE, Class=TRUE){
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
			cptsSig <- data.frame('cpts'=cpts[1],'T'=T,'thresholdValue'=thresholdValue,'significant'=TRUE)
		}else{
			cpts <- c(n)
			cptsSig <- data.frame('cpts'=which.max(testStat),'T'=T,'thresholdValue'=thresholdValue,'significant'=FALSE)
		}
	}else{
		if(is.numeric(numCpts)){
			thresholdValue <- 0
			cptsSig <- binSeg(X,method='RATIO',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=numCpts)
		}else{
			cptsSig <- binSeg(X,method='RATIO',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1)
		}
		cpts <- c(cptsSig$cpts[cptsSig$significant],n)
	}
	if(Class==FALSE){
		return(cpts)
	}else{
		return(classInput(X=X,cpts=cpts,method='Ratio',numCpts=numCpts,cptsSig=cptsSig,threshold=threshold,thresholdValue=unique(cptsSig$thresholdValue),msl=msl))
	}
}

