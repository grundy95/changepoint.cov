#' @include cptCovClass.R
NULL

#' Extension of S4 cptCov class for the Subspace method 
#'
#' This S4 class extends the \code{\link{cptCov}} for detailed analysis of results from the Subspace method. NOTE only methods and slots that aren't included in the \code{\link{cptCov}} are described here. Your first port of call should be the \code{\link{cptCov}} documentation.
#' 
#' @slot q Assumed dimension of latent subspace
#' @slot nperm Number of permutations used in the permutation test if this was the chosen method for calculating the threshold.
setClass("cptCovSubspace",slots=list(q='numeric',nperm='numeric'),prototype=list(method='Subspace',version=as(packageVersion("changepoint.cov"),'character'),date=date(),method=NULL),contains='cptCov')


#' @describeIn cptCovSubspace Shows the cptCovSubspace object
#'
#' Shows the cptCov object and all the slots it contains. Also displays the summary of the object.
setMethod("show","cptCovSubspace",function(object){
		  cat('Class, cptCovSubspace     : Covariance Changepoint Object\n')
		  cat('		    	         : S4 class containing ',length(attributes(object))-1,' slots with names\n')
		  cat('	                          ',names(attributes(object))[1:(length(attributes(object))-1)],'\n\n')
		  cat('Created on    	          : ',object@date,'\n\n')
		  cat('Summary(.)    	          :\n---------------\n')
		  summary(object)
})
