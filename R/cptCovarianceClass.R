#' cptCovariance: an \code{S4} class for a covariance changepoint object
#'
#' Contains data and information required for further changepoint analysis,
#' summaries and plotting.
#'
#' @slot dataset An n by p matrix of the data.
#' @slot cpts A numeric vector containing the identified changepoints.
#' @slot method Character containing covariance changepoint method used.
#' @slot numCpts Either 'AMOC' for at most one changepoint; 'BinSeg' for a binary segmentation approach to detect multiple changepoints; or a positive integer specifying the number of changepoints.
#' @slot cptsSig Data frame containing the changepoint locations along with their associated test statistic; threshold; and whether or not they were deemed significant.
#' @slot threshold Character containing the method used for generating the threshold.
#' @slot thresholdValue Threshold value used to determine significant changepoints. If permutation test within subspace method is used then a vector of threshold values is returned.
#' @slot msl Minimum segment length between changepoints.
#' @slot subspaceDim Assumed subspace dimension. Only used for Subspace method.
#' @slot nperm Numeric value of number of permutations used in permutation test. Only used for subspace method when threshold is "PermTest".
#' @slot LRCov Character describing the long-run covariance estimator used or a matrix containing the long-run covariance estimate. Only used for CUSUM method.
#' @slot covEst List containing the sample covariance for each segment.
#' @slot subspaceEst List containing a basis of the subspace for each segment. Only used for Subspace method.
#' @slot date Creation date of the object.
#' @slot version Version of the cpt.covariance package used.
#'
#' @examples
#' # Create new cptCovaraince object
#' out <- new('cptCovariance',
#'      dataset=matrix(rnorm(300),ncol=3),
#'			cpts=c(50,100),
#'			method='Ratio',
#'			numCpts='AMOC',
#'			cptsSig=data.frame('cpts'=50,'T'=33.3,'thresholdValue'=30,significant=TRUE),
#'			threshold='Manual',
#'			thresholdValue=30,
#'			msl=20)
#'
#' # Summarize cptCovariance object
#' summary(out)
#'
#' # Show cptCovariance object
#' show(out)
#'
#' # Plot cptCovariance object
#' plot(out)
#'
#' # Show significant changepoints
#' cptsSig(out)
#'
#' # Estimate covariance matrices in each segment
#' covEst(out)
#'
#' @import methods
#' @export
setClass("cptCovariance",slots=list(dataset='matrix',cpts='numeric',method='character',msl='numeric',numCpts='ANY',threshold='character',thresholdValue='numeric',cptsSig='data.frame',subspaceDim='numeric',nperm='numeric',LRCov='ANY',covEst='list',subspaceEst='list',date='character',version='character'),prototype=list(subspaceDim=NA_real_,nperm=NA_real_,LRCov=NA_character_,covEst=list(NA_real_),subspaceEst=list(NA_real_),version=as(packageVersion("changepoint.cov"),'character'),date=date(),method=NULL))

#' @describeIn cptCovariance Summarises the cptCovariance object
#'
#' @param object An object of S4 class \code{\linkS4class{cptCovariance}}
#' @export
setMethod("summary","cptCovariance",function(object){
		  cat('Created using changepoint.cov version',object@version,'\n')
		  cat('Method			    : ',object@method,'\n')
		  cat('Multiple changepoint method : ',numCpts(object),'\n')
		  cat('Minimum segment length      : ',object@msl,'\n')
		  cat('Changepoints		    : ',object@cpts,'\n')
})

#' @describeIn cptCovariance Shows the cptCovariance object
#'
#' @param object An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @export
setMethod("show","cptCovariance",function(object){
		  cat('Class, cptCovariance  : Covariance Changepoint x\n')
		  cat('		      : S4 class containing ',length(attributes(object))-1,' slots with names\n')
		  cat('	               ',names(attributes(object))[1:(length(attributes(object))-1)],'\n\n')
		  cat('Created on      : ',object@date,'\n\n')
		  cat('Summary(.)      :\n---------------\n')
		  summary(object)
})


#' @describeIn cptCovariance Plotting method for cptCovariance object. Returns a \code{\link[ggplot2]{ggplot}} object which can be manipulated as required
#'
#' @param object An object of S4 class \code{\linkS4class{cptCovariance}}
#' @param x x
#' @param y y
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod("plot", "cptCovariance", function(object, x=c(), y=c()){
		  p <- ncol(object@dataset)
		  covs <- covEst(object)

		  segCovs <- data.frame('Segment'=1,'p1'=rep(1:p,each=p),'p2'=rep(1:p,p),'Value'=as.vector(covs[[1]]))
		  if(length(cpts(object))>1){
			  for(i in 2:length(cpts(object))){
				  segCovs <- rbind(segCovs,data.frame('Segment'=i,'p1'=rep(1:p,each=p),'p2'=rep(1:p,p),'Value'=as.vector(covs[[i]])))
			  }
		  }
		  covPlot <- ggplot(segCovs,aes(x=.data$p1,y=.data$p2,fill=.data$Value))+
			  geom_tile()+
			  facet_grid(.~Segment,labeller='label_both')+
			  scale_fill_viridis(name="Value")+
			  xlab("p")+ylab("p")
		  return(covPlot)
})




setGeneric("covEst", function(object){
  standardGeneric("covEst")
})

#' @describeIn cptCovariance Returns covariance estimates for each segment
#'
#' @param object An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases covEst
#'
#' @export
setMethod("covEst","cptCovariance",function(object){
		  X <- object@dataset
		  covs <- list()
		  cpts <- c(0,cpts(object))
		  for(i in 1:(length(cpts)-1)){
			  covs[[i]] <- cov(X[(cpts[i]+1):cpts[i+1],])
		  }
		  return(covs)
})

setGeneric("subspaceEst", function(object){
  standardGeneric("subspaceEst")
})

#' @describeIn cptCovariance Returns a basis of the subspace estimates for each segment
#'
#' @param object An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases subspaceEst
#'
#' @export
setMethod("subspaceEst","cptCovariance",function(object){
		  if(method(object)!='Subspace'){
			  stop("Subspace estimation only possible for method='Subspace'")
		  }else{
			  X <- object@dataset
			  q <- subspaceDim(object)
			  subspace <- list()
			  cpts <- c(0,cpts(object))
			  for(i in 1:(length(cpts)-1)){
				 subspace[[i]] <- eigen(cov(X[(cpts[i]+1):cpts[i+1],]),symmetric=TRUE)$vectors[,1:q]
			  }
		  }
		  return(subspace)
})




#' Retrieval Functions - Generic
#'
#' @param x object of class \code{\linkS4class{cptCovariance}}
#' @name retrievalGeneric
NULL

#' Replacement Functions - Generic
#'
#' @param x object of class \code{\linkS4class{cptCovariance}}
#' @param value value
#' @name replacementGeneric
NULL

#' Retrieval Functions - Method
#'
#' @param x object of class \code{\linkS4class{cptCovariance}}
#' @name retrievalMethod
NULL

#' Replacement Functions - Method
#'
#' @param x object of class \code{\linkS4class{cptCovariance}}
#' @param value value
#' @name replacementMethod
NULL

#' @rdname retrievalGeneric
#' @export
setGeneric("dataset", function(x){
  standardGeneric("dataset")
})

#' @rdname replacementGeneric
#' @export
setGeneric("dataset<-", function(x, value){
  standardGeneric("dataset<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("dataset", "cptCovariance", function(x){
  x@dataset
})

#' @rdname replacementMethod
#' @export
setMethod("dataset<-", "cptCovariance", function(x, value){
  x@dataset <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("cptsSig", function(x){
  standardGeneric("cptsSig")
})

#' @rdname replacementGeneric
#' @export
setGeneric("cptsSig<-", function(x, value){
  standardGeneric("cptsSig<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("cptsSig", "cptCovariance", function(x){
  x@cptsSig
})

#' @rdname replacementMethod
#' @export
setMethod("cptsSig<-", "cptCovariance", function(x, value){
  x@cptsSig <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("cpts", function(x){
  standardGeneric("cpts")
})

#' @rdname replacementGeneric
#' @export
setGeneric("cpts<-", function(x, value){
  standardGeneric("cpts<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("cpts", "cptCovariance", function(x){
  x@cpts
})

#' @rdname replacementMethod
#' @export
setMethod("cpts<-", "cptCovariance", function(x, value){
  x@cpts <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("method", function(x){
  standardGeneric("method")
})

#' @rdname replacementGeneric
#' @export
setGeneric("method<-", function(x, value){
  standardGeneric("method<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("method", "cptCovariance", function(x){
  x@method
})

#' @rdname replacementMethod
#' @export
setMethod("method<-", "cptCovariance", function(x, value){
  x@method <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("msl", function(x){
  standardGeneric("msl")
})

#' @rdname replacementGeneric
#' @export
setGeneric("msl<-", function(x, value){
  standardGeneric("msl<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("msl", "cptCovariance", function(x){
  x@msl
})

#' @rdname replacementMethod
#' @export
setMethod("msl<-", "cptCovariance", function(x, value){
  x@msl <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("numCpts", function(x){
  standardGeneric("numCpts")
})

#' @rdname replacementGeneric
#' @export
setGeneric("numCpts<-", function(x, value){
  standardGeneric("numCpts<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("numCpts", "cptCovariance", function(x){
  if(is.numeric(x@numCpts)){
    return(paste0("Manual - ",x@numCpts," changepoints"))
  }else{
    return(x@numCpts)
  }
})

#' @rdname replacementMethod
#' @export
setMethod("numCpts<-", "cptCovariance", function(x, value){
  x@numCpts <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("threshold", function(x){
  standardGeneric("threshold")
})

#' @rdname replacementGeneric
#' @export
setGeneric("threshold<-", function(x, value){
  standardGeneric("threshold<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("threshold", "cptCovariance", function(x){
  x@threshold
})

#' @rdname replacementMethod
#' @export
setMethod("threshold<-", "cptCovariance", function(x, value){
  x@threshold <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("thresholdValue", function(x){
  standardGeneric("thresholdValue")
})

#' @rdname replacementGeneric
#' @export
setGeneric("thresholdValue<-", function(x, value){
  standardGeneric("thresholdValue<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("thresholdValue", "cptCovariance", function(x){
  x@thresholdValue
})

#' @rdname replacementMethod
#' @export
setMethod("thresholdValue<-", "cptCovariance", function(x, value){
  x@thresholdValue <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("subspaceDim", function(x){
  standardGeneric("subspaceDim")
})

#' @rdname replacementGeneric
#' @export
setGeneric("subspaceDim<-", function(x, value){
  standardGeneric("subspaceDim<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("subspaceDim", "cptCovariance", function(x){
  if(toupper(x@method)!='SUBSPACE'){
    stop("subspaceDim is only a valid slot for method='Subspace'")
  }else{
    return(x@subspaceDim)
  }
})

#' @rdname replacementMethod
#' @export
setMethod("subspaceDim<-", "cptCovariance", function(x, value){
  x@subspaceDim <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("nperm", function(x){
  standardGeneric("nperm")
})

#' @rdname replacementGeneric
#' @export
setGeneric("nperm<-", function(x, value){
  standardGeneric("nperm<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("nperm", "cptCovariance", function(x){
  if(!((toupper(x@method)=='SUBSPACE')&&(toupper(x@threshold)=='PERMTEST'))){
    stop("nperm is only a valid slot when using the permutation test within method='Subspace'")
  }else{
    return(x@nperm)
  }
})

#' @rdname replacementMethod
#' @export
setMethod("nperm<-", "cptCovariance", function(x, value){
  x@nperm <- value
})

#' @rdname retrievalGeneric
#' @export
setGeneric("LRCov", function(x){
  standardGeneric("LRCov")
})

#' @rdname replacementGeneric
#' @export
setGeneric("LRCov<-", function(x, value){
  standardGeneric("LRCov<-")
})

#' @rdname retrievalMethod
#' @export
setMethod("LRCov", "cptCovariance", function(x){
  if(toupper(x@method)!='CUSUM'){
    stop("LRCov is only a valid slot for method='CUSUM'")
  }else{
    return(x@LRCov)
  }
})

#' @rdname replacementMethod
#' @export
setMethod("LRCov<-", "cptCovariance", function(x, value){
  x@LRCov <- value
})

