#' Inputs results into S4 class
#'
#' Function for inputting results into the S4 class \code{\link{cptCov}}. 
#'
classInput <- function(X,cpts,method,noCpts,testStat,threshold,thresholdValue,msl,q=0,nperm=0,LRCov='NA',statType='NA'){
	if(method=='Subspace'){
		ans <- new('cptCov',data=X,
                           cpts=cpts,
			   method=method,
	                   noCpts=noCpts,
	         	   testStat=testStat,
			   threshold=threshold,
			   thresholdValue=thresholdValue,
			   msl=msl,
			   q=q,
			   nperm=nperm)
	}else if(method=='CUSUM'){
		ans <- new('cptCov',data=X,
			   cpts=cpts,
			   method=method,
			   noCpts=noCpts,
			   testStat=testStat,
			   threshold=threshold,
			   thresholdValue=thresholdValue,
			   msl=msl,
			   LRCov=LRCov,
			   statType=statType)
	}else if(method=='Fisher'){
		ans <- new('cptCov',data=X,
			   cpts=cpts,
			   method=method,
			   noCpts=noCpts,
			   testStat=testStat,
			   threshold=threshold,
			   thresholdValue=thresholdValue,
			   msl=msl)
	}
	return(ans)
}
