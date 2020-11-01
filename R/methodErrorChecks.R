#' Error checking for Ratio Method
#'
#' DEVELOPER USE ONLY. This function checks the user inputs to make sure they are all valid.
#'
#' @inheritParams cptRatio
#'
#' @seealso \code{\link{cptRatio}}
#'
#' @keywords internal
ratioErrorChecks <- function(X,threshold,numCpts,thresholdValue,msl,Class){
	dataErrorChecks(X)
	thresholdErrorChecks(threshold,thresholdValue,method='cptCov')
	numCptsErrorChecks(numCpts,method="cptCov")
	mslErrorChecks(X,msl)
	classErrorChecks(Class)

}

#' Error checking for CUSUM Method
#'
#' DEVELOPER USE ONLY. This function checks the user inputs to make sure they are all valid
#'
#' @inheritParams cptCUSUM
#'
#' @seealso \code{\link{cptCUSUM}}
#' @keywords internal
cusumErrorChecks <- function(X,threshold,numCpts,LRCov,thresholdValue,msl,Class){
	dataErrorChecks(X)
	thresholdErrorChecks(threshold,thresholdValue,method='cptCov')
	numCptsErrorChecks(numCpts,method='cptCov')
	mslErrorChecks(X,msl)
	classErrorChecks(Class)

	#LRCov checks
	p <- ncol(X)
	delta <- (p*(p+1))/2
	if(!(is.character(LRCov)||is.matrix(LRCov))){
		stop("LRCov not identified: see ?cptCov for valid entries to LRCov")
	}
	if((is.matrix(LRCov))&&(!((ncol(LRCov)==delta)&&nrow(LRCov==delta)))){
		stop("Dimension of manual LRCov is not compatible with data")
	}	
	if(is.character(LRCov)){
		LRCov <- toupper(LRCov)
		if((LRCov!="BARTLETT")&&(LRCov!="EMPIRICAL")){
			stop("LRCov not identified: see ?cptCov for valid entries to LRCov")
		}
	}
}

#' Error checking for Subspace Method
#'
#' DEVELOPER USE ONLY. This function checks the user inputs to make sure they are all valid
#'
#' @inheritParams cptSubspace
#'
#' @seealso \code{\link{cptSubspace}}
#'
#' @keywords internal
subspaceErrorChecks <- function(X,subspaceDim,threshold,numCpts,thresholdValue,msl,nperm,Class){
	dataErrorChecks(X)
	thresholdErrorChecks(threshold,thresholdValue,method='cptSubspace')
	numCptsErrorChecks(numCpts,method="cptSubspace")
	mslErrorChecks(X,msl)
	classErrorChecks(Class)

	#subspaceDim checks
	if(!is.numeric(subspaceDim)){
		stop("Subspace dimension should be a single positive integer")
	}
	if((length(subspaceDim)!=1)||(subspaceDim%%1!=0)||(subspaceDim<1)){
		stop("Subspace dimension should be a single positive integer")
	}
	if(subspaceDim>ncol(X)){
		stop("Subspace dimension, subspaceDim, should be smaller than time series dimension, p.")
	}


	#nperm checks
	if(!is.numeric(nperm)){
		stop("Number of permutations should be a single positive integer")
	}
	if((toupper(threshold)!="PERMTEST")&&(nperm!=200)){
		warning("nperm is only used with threshold='PermTest'")
	}
	if((length(nperm)!=1)||(nperm%%1!=0)||(nperm<1)){
		stop("Number of permutations should be a single positive integer")
	}
	if(nperm<50){
		warning("Number of permutations is small - threshold may be unreliable")
	}
	if(nperm>(nrow(X)*5)){
		warning("Number of permutations is large - method may take substantial time to run")
	}
}





	





