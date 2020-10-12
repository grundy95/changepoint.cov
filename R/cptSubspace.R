#' Detecting changes in subspace
#'
#' This method allows the detection of subspace changes in multivariate time series data.
#'
#' @param X Data matrix of dimension n by p
#' @param q Dimension of the latent subspace
#' @param threshold Threshold choice for determining significance of changepoints. Choices include:
#' \itemize{
#'	\item "PermTest" - Permutation test is performed using the number of permutations and significance level contained in the nperm and thresholdValue parameters respectively
#'	\item "Manual" - A user chosen threshold is used which is contained in the thresholdValue argument.
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
#' @param Class Logical. If TRUE then an S4 class is returned. Else just the estimated changepoints are returned.
#'
#' @examples
#' set.seed(1)
#' dataAMOC <- subspaceDataGeneration(n=100,p=20,q=5,tau=50,changeSize=0.5*sqrt(5))$data
#' dataMultipleCpts <- subspaceDataGeneration(n=200,p=20,q=5,tau=c(50,100,150),
#'						changeSize=0.4*sqrt(5))$data
#' 
#' set.seed(1)
#' ansSubspace <- cptSubspace(X=dataAMOC,q=5,nperm=100)
#' summary(ansSubspace)
#' 
#' ansSubspace2 <- cptSubspace(X=dataMultipleCpts,q=5,threshold='Manual',numCpts='BinSeg',
#'				thresholdValue=30,msl=30)
#' summary(ansSubspace2)
#' 
#' ansSubspace3 <- cptSubspace(X=dataMultipleCpts,q=5,numCpts=3,msl=30)
#' summary(ansSubspace3)
#'
#' @return An object of S4 class \code{\link{cptCovariance-class}} is returned. If Class="FALSE", just the vector of changepoints are returned.
#' @export
cptSubspace <- function(X,q,threshold='PermTest',numCpts='AMOC',thresholdValue=0.05,msl=ncol(X),nperm=200,Class=TRUE){
	if(is.data.frame(X)){
		X <- as.matrix(X)
	}
	subspaceErrorChecks(X=X,q=q,threshold=threshold,numCpts=numCpts,thresholdValue=thresholdValue,msl=msl,nperm=nperm,Class=Class)
	threshold <- toupper(threshold)
	if(is.character(numCpts)){
		numCpts <- toupper(numCpts)
	}
	n <- nrow(X)
	p <- ncol(X)
	if(numCpts=='AMOC'){
		testStat <- subspaceTestStat(X,q,msl)
		T <- max(testStat,na.rm=TRUE)
		if(threshold=='PERMTEST'){
			thresholdValue <- permutationTest(X,q,msl,thresholdValue,nperm)
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
			cptsSig <- binSeg(X,q=q,method='SUBSPACE',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=numCpts)
		}else{
			if(threshold=='PERMTEST'){
				alpha <- thresholdValue
				cptsSig <- binSeg(X,method='SUBSPACE',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1,q=q,alpha=alpha,nperm=nperm)
			}else{
				cptsSig <- binSeg(X,method='SUBSPACE',msl=msl,threshold=threshold,thresholdValue=thresholdValue,m=-1,q=q)
			}
		}
		cpts <- c(cptsSig$cpts[cptsSig$significant],n)
	}
	if(Class==FALSE){
		return(cpts)
	}else{
		return(classInput(X=X,cpts=cpts,method='Subspace',numCpts=numCpts,testStat=cptsSig$T,threshold=threshold,thresholdValue=cptsSig$thresholdValue,msl=msl,q=q,nperm=nperm))
	}
}

