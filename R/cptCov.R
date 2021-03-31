#' Covariance changepoint detection
#'
#' Finds covariance changepoints within multivariate time series data using
#' either the Ratio method \insertCite{Ryan2020}{changepoint.cov} or a
#' CUSUM based method \insertCite{Aue2009}{changepoint.cov}.
#'
#' This is a simple wrapper function for the functions
#' \code{\link{cptRatio}} and \code{\link{cptCUSUM}}. If no method is specified
#' then the method used will depend on the dimension of the time series. For
#' p<10, the CUSUM method will be used and for p>=10 the Ratio method will be
#' used.
#'
#' @param X Data matrix of dimension n by p.
#' @param method Covariance changepoint method to be used. Choice of "Ratio" or "CUSUM".
#' @param threshold Threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "Asymptotic" - Uses the asymptotic threshold derived for each method. For Ratio method the threshold is log(n). For CUSUM method the threshold is the specified quantile of the standard Normal distribution. The quantile is set by the argument thresholdValue.
#'	\item "Manual"- A user chosen threshold which is contained in the thresholdValue argument. NOTE the normalized test statistics will be compared to the set thresholds - see details for more information.
#'}
#' If numCpts is numeric then the threshold is not used as the number of changepoints is known.
#' @param numCpts Number of changepoints in the data. Choices include:
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item Numeric - User specified number of changepoints.
#' }
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param LRCov The long-run covariance estimator to be used for CUSUM method. Currently, only "Bartlett" and "Empirical" are supported. Alternatively, a matrix containing the long-run covariance estimate can be inputted.
#' @param thresholdValue Either the manual threshold value when threshold="Manual" or the (1-thresholdValue)-quantile of asymptotic distribution of the test statistic when threshold="Asymptotic".
#' @param Class Logical. If TRUE then an S4 class is returned. If FALSE the estimated changepoints are returned.
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#'
#' @references
#' \insertRef{Ryan2020}{changepoint.cov}
#'
#' \insertRef{Aue2009}{changepoint.cov}
#'
#' @seealso \code{\link{cptRatio}}, \code{\link{cptCUSUM}}, \code{\linkS4class{cptCovariance}}, \code{\link{wishartDataGeneration}}
#'
#' @examples
#' set.seed(1)
#' dataAMOC <- wishartDataGeneration(n=100,p=3,tau=50)$data
#' dataMultipleCpts <- wishartDataGeneration(n=200,p=3,tau=c(50,100,150))$data
#'
#' ansRatio <- cptCov(X=dataAMOC,method="Ratio")
#' summary(ansRatio)
#' plot(ansRatio)
#'
#' ansCUSUM <- cptCov(X=dataAMOC,method='CUSUM')
#' show(ansCUSUM)
#'
#' ansRatio2 <- cptCov(X=dataMultipleCpts,method='Ratio',threshold='Manual',numCpts='BinSeg',
#'			msl=10,thresholdValue=20)
#' summary(ansRatio2)
#' cptsSig(ansRatio2)
#' plot(ansRatio2)
#'
#' ansCUSUM2 <- cptCov(X=dataAMOC,method='CUSUM',numCpts=3,
#'			msl=15,LRCov='Empirical')
#' summary(ansCUSUM2)
#' cptsSig(ansCUSUM2)
#'
#' @include cptCovariance-class.R
#'
#' @export
cptCov <- function(X,method=c("Ratio","CUSUM"),threshold="Asymptotic",numCpts='AMOC',msl=2*ncol(X),thresholdValue=0.05,LRCov='Bartlett',Class=TRUE){
	if(!is.character(method)){
		stop("method not recognized: Please choose between 'Ratio' and 'CUSUM'")
	}
	if((length(method)>1)&&(method[1]=="Ratio")&&(method[2]=="CUSUM")){
		if(ncol(X)>10){
			method <- "Ratio"
			warning("no method was chosen. As p>10 the Ratio method will be implemented")
		}else{
			method <- "CUSUM"
			warning("no method was chosen. As p<=10 the CUSUM method will be implemented")
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
		if(is.character(numCpts)){
			numCpts <- toupper(numCpts)
		}
		ans <- cptRatio(X=X,threshold=threshold,numCpts=numCpts,msl=msl,thresholdValue=thresholdValue,errorCheck=FALSE,Class=Class)
	}else{
		cusumErrorChecks(X=X,threshold=threshold,numCpts=numCpts,msl=msl,LRCov=LRCov,thresholdValue=thresholdValue,Class=Class)
		threshold <- toupper(threshold)
		if(is.character(numCpts)){
			numCpts <- toupper(numCpts)
		}
		if(is.character(LRCov)){
			LRCov <- toupper(LRCov)
		}
		ans <- cptCUSUM(X=X,threshold=threshold,numCpts=numCpts,msl=msl,LRCov=LRCov,thresholdValue=thresholdValue,errorCheck=FALSE,Class=Class)
	}
	return(ans)
}



