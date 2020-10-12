#' Inputs results into S4 class
#'
#' Function for inputting results into the S4 class \code{\link{cptCovariance-class}}. 
#'
#' @inheritParams cptCov
#' @inheritParams cptSubspace
#' @param cptsSig Data frame containing the changepoint locations along with their associated test statistics; threshold; and whether or not they were deemed significant  
#' @param cpts Vector of changepoint locations
#'
#' @return S4 class of type \code{\link{cptCovariance-class}}
#'
#' @importFrom methods new
classInput <- function(X,cpts,method,numCpts,cptsSig,threshold,thresholdValue,msl,q=0,nperm=0,LRCov='NA'){
	if(method=='Subspace'){
		ans <- new('cptCovariance',data=X,
                           cpts=cpts,
			   method=method,
	                   numCpts=numCpts,
	         	   cptsSig=cptsSig,
			   threshold=threshold,
			   thresholdValue=thresholdValue,
			   msl=as.integer(msl),
			   q=q,
			   nperm=nperm)
	}else if(method=='CUSUM'){
		ans <- new('cptCovariance',data=X,
			   cpts=cpts,
			   method=method,
			   numCpts=numCpts,
			   cptsSig=cptsSig,
			   threshold=threshold,
			   thresholdValue=thresholdValue,
			   msl=as.integer(msl),
			   LRCov=LRCov)
	}else if(method=='Ratio'){
		ans <- new('cptCovariance',data=X,
			   cpts=cpts,
			   method=method,
			   numCpts=numCpts,
			   cptsSig=cptsSig,
			   threshold=threshold,
			   thresholdValue=thresholdValue,
			   msl=as.integer(msl))
	}
	return(ans)
}
