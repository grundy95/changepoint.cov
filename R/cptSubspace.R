#' Detecting changes in subspace
#'
#' Implements the \insertCite{Grundy2020;textual}{changepoint.cov} method for detecting changes in subspace in multivariate time series data. This method is aimed at time series that lie in a lower-dimensional latent subspace.
#'
#' Subspace changepoint detection is aimed at time series where we assume the data lies in a low-dimensional subspace; meaning there are a q dominating eigenvalues in the covariance matrix, where q is the assumed subspace dimension. This function calculates the test statistic, $T$ described in \insertCite{Grundy2020;textual}{changepoint.cov}. A data driven threshold is recommended by using the permutation test to determine the significance of changepoints. Note that this is a data driven threshold and will therefore vary in each calculation. The calculation of the threshold via the permutation test can also be computationally expensive especially for long time series (n>1000) or a large number of permutations (nperm>1000). The number of permutations should be altered to reflect the length of the time series - the longer the time series the more permutations may be necessary. If multiple changepoints are possible then the Binary Segmentation is implemented however only one segment will be tested at each iteration in order to control the type 1 error. This method is only recommended if the data is assumed to lie in a low-dimension subspace and the dimensionality of this subspace is known. If this is not the case we recommend using the \code{\link{cptCov}} function and one of the contained methods within.
#'
#' @param X Data matrix of dimension n by p.
#' @param subspaceDim Dimension of the latent subspace.
#' @param threshold Threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "PermTest" - Permutation test is performed using the number of permutations and significance level contained in the nperm and thresholdValue parameters respectively.
#'	\item "Manual" - A user chosen threshold is used, which is contained in the thresholdValue argument.
#' }
#' If numCpts is numeric then the threshold is not used as the number of changepoints is known.
#' @param numCpts Number of changepoints in the data. Choices include: 
#' \itemize{
#' 	\item "AMOC" - At Most One Changepoint; test to see if the data contains a single changepoint or not.
#'	\item "BinSeg"- Binary segmentation is performed to detect multiple changepoints.
#'	\item Numeric - User specified number of changepoints.
#' }
#' @param thresholdValue Either the significance level of the permutation test when using threshold="PermTest" or the user defined threshold when using threshold="Manual".
#' @param msl Minimum segment length allowed between the changepoints. NOTE this should be greater than or equal to p, the dimension of the time series.
#' @param nperm Only required for threshold="PermTest". Number of permutations to use in the permutation test.
#' @param Class Logical. If TRUE then an S4 class is returned. Else the estimated changepoints are returned.
#'
#' @references 
#' \insertRef{Grundy2020}{changepoint.cov}
#'
#' @seealso \code{\link{cptCov}}, \code{\linkS4class{cptCovariance}}, \code{\link{subspaceDataGeneration}}, \code{\link{permutationTest}}
#'
#'
#' @examples
#' set.seed(1)
#' dataAMOC <- subspaceDataGeneration(n=100,p=20,subspaceDim=5,tau=50,changeSize=0.5*sqrt(5))$data
#' dataMultipleCpts <- subspaceDataGeneration(n=200,p=20,subspaceDim=5,tau=c(50,100,150),
#'						changeSize=0.4*sqrt(5))$data
#' 
#' set.seed(1)
#' ansSubspace <- cptSubspace(X=dataAMOC,subspaceDim=5,nperm=100)
#' summary(ansSubspace)
#' subspaceEst(ansSubspace)
#' 
#' ansSubspace2 <- cptSubspace(X=dataMultipleCpts,subspaceDim=5,threshold='Manual',numCpts='BinSeg',
#'				thresholdValue=30,msl=30)
#' summary(ansSubspace2)
#' cptsSig(ansSubspace2)
#' 
#' ansSubspace3 <- cptSubspace(X=dataMultipleCpts,subspaceDim=5,numCpts=3,msl=30)
#' summary(ansSubspace3)
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", the vector of changepoints are returned.
#'
#' @include cptCovarianceClass.R
#'
#' @export
cptSubspace <- function(X,subspaceDim,threshold='PermTest',numCpts='AMOC',thresholdValue=0.05,msl=ncol(X),nperm=200,Class=TRUE){
	if(is.data.frame(X)){
		X <- as.matrix(X)
	}
	subspaceErrorChecks(X=X,subspaceDim=subspaceDim,threshold=threshold,numCpts=numCpts,thresholdValue=thresholdValue,msl=msl,nperm=nperm,Class=Class)
	threshold <- toupper(threshold)
	if(is.character(numCpts)){
		numCpts <- toupper(numCpts)
	}
	n <- nrow(X)
	p <- ncol(X)
	if(numCpts=='AMOC'){
		testStat <- subspaceTestStat(X,subspaceDim,msl)
		T <- max(testStat,na.rm=TRUE)
		if(threshold=='PERMTEST'){
			thresholdValue <- permutationTest(X,subspaceDim,msl,thresholdValue,nperm)
		}else if(threshold=='MANUAL'){
			nperm <- 0
		}
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
			cptsSig <- binSeg(X,subspaceDim=subspaceDim,method='SUBSPACE',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=numCpts)
		}else{
			if(threshold=='PERMTEST'){
				alpha <- thresholdValue
				cptsSig <- binSeg(X,method='SUBSPACE',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1,subspaceDim=subspaceDim,alpha=alpha,nperm=nperm)
			}else{
				cptsSig <- binSeg(X,method='SUBSPACE',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1,subspaceDim=subspaceDim)
			}
		}
		cpts <- c(cptsSig$cpts[cptsSig$significant],n)
	}
	if(Class==FALSE){
		return(cpts)
	}else{
		return(classInput(X=X,cpts=cpts,method='Subspace',numCpts=numCpts,cptsSig=cptsSig,threshold=threshold,thresholdValue=unique(cptsSig$thresholdValue),msl=msl,subspaceDim=subspaceDim,nperm=nperm))
	}
}

