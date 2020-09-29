#' An S4 class for a covariance changepoint object
#' 
#' @slot data An n by p matrix of the data
#' @slot cpts A numeric vector containing the identified changepoints
#' @slot method Character containing covariance changepoint method used
#' @slot noCpts Either 'AMOC' for at most one changepoint; 'BinSeg' for a binary segmentation approach to dect multiple changepoints; or a positive integer specifying the number of changes
#' @slot testStat Numeric containing the test statistic for each changepoint
#' @slot threshold Character containing the method used for generating the threshold
#' @slot thresholdValue Numeric value of threshold
#' @slot msl Numeric containing the minimum segment length between changepoints
#' @slot q Numeric value of subspace dimension. Only used for Subspace method
#' @slot nperm Numeric value of number of permutations used in permutation test. Only used for subspace method when threshold is "PermTest"
#' @slot LRCov Character describing the long-run covariance estimator used. Only used for Aue method.
#' @slot statType Character describing which of the two test statistics were used. Only used for Aue method.
#' @slot date Creation date of the object
#' @slot version Version of the cpt.covariance used
#' @export
setClass("cptCov",slots=list(data='matrix',cpts='numeric',method='character',msl='numeric',noCpts='ANY',threshold='character',thresholdValue='numeric',testStat='numeric',q='numeric',nperm='numeric',LRCov='character',statType='character',date='character',version='character'),prototype=list(q=0,nperm=0,LRCov='NA',statType='NA',version=as(packageVersion("changepoint.cov"),'character'),date=date(),method=NULL))

#' @describeIn cptCov Summarises the cptCov object
#'
#' Summaries the results of the covariance changepoint analysis contained in the cptCov object.
setMethod("summary","cptCov",function(object){
		  cat('Created using changepoint.cov version',object@version,'\n')
		  cat('Method			    : ',object@method,'\n')
		  cat('Multiple changepoint method : ',object@noCpts,'\n')
		  cat('Minimum segment length      : ',object@msl,'\n')
		  cat('Changepoints		    : ',object@cpts,'\n')
})

#' @describeIn cptCov Shows the cptCov object
#'
#' Shows the cptCov object and all the slots it contains. Also displays the summary of the object.
setMethod("show","cptCov",function(object){
		  cat('Class, cptCov  : Covariance Changepoint Object\n')
		  cat('		      : S4 class containing ',length(attributes(object))-1,' slots with names\n')
		  cat('	               ',names(attributes(object))[1:(length(attributes(object))-1)],'\n\n')
		  cat('Created on      : ',object@date,'\n\n')
		  cat('Summary(.)      :\n---------------\n')
		  summary(object)
})
