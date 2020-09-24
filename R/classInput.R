#' @include cptCovClass.R cptCovSubspaceClass.R
NULL

#' Inputs results into S4 class
#'
#' Function for inputting results into the S4 class \code{\link{cptCov}}. 
#'

classInput <- function(X,cpts,method,noCpts,testStat,threshold,thresholdValue,msl,q,nperm){
	if(method=='Subspace'){
		ans <- new('cptCovSubspace',data=X,
				      cpts=cpts,
				      method=method,
				      noCpts=noCpts,
				      testStat=testStat,
				      threshold=threshold,
				      thresholdValue=thresholdValue,
				      msl=msl,
				      q=q,
				      nperm=nperm)
	}
	return(ans)
}
