#' An S4 class for a covariance changepoint object
#' 
#' @slot data An n by p matrix of the data
#' @slot cpts A numeric vector containing the identified changepoints
#' @slot method Character containing covariance changepoint method used
#' @slot numCpts Either 'AMOC' for at most one changepoint; 'BinSeg' for a binary segmentation approach to dect multiple changepoints; or a positive integer specifying the number of changes
#' @slot cptsSig Data frame containing the changepoint locations along with their associated test statistic; threshold; and whether or not they were deemed significant
#' @slot threshold Character containing the method used for generating the threshold
#' @slot thresholdValue Threshold value used to determine significant changepoints. If permutation test within subspace method is used then a vector of threshold values is returned. 
#' @slot msl Numeric containing the minimum segment length between changepoints
#' @slot subspaceDim Numeric value of subspace dimension. Only used for Subspace method
#' @slot nperm Numeric value of number of permutations used in permutation test. Only used for subspace method when threshold is "PermTest"
#' @slot LRCov Character describing the long-run covariance estimator used. Only used for Aue method.
#' @slot date Creation date of the object
#' @slot version Version of the cpt.covariance used
#'
#' @examples
#' ans <- new('cptCovariance',data=matrix(rnorm(300),ncol=3),
#'			cpts=c(50,100),
#'			method='Ratio',
#'			numCpts='AMOC',
#'			cptsSig=data.frame('cpts'=50,'T'=33.3,'thresholdValue'=30,significant=TRUE),
#'			threshold='Manual',
#'			thresholdValue=30,
#'			msl=20)
#' summary(ans)
#' show(ans)
#'
#' @export
setClass("cptCovariance",slots=list(data='matrix',cpts='numeric',method='character',msl='numeric',numCpts='ANY',threshold='character',thresholdValue='numeric',cptsSig='data.frame',subspaceDim='numeric',nperm='numeric',LRCov='character',date='character',version='character'),prototype=list(subspaceDim=NA_real_,nperm=NA_real_,LRCov=NA_character_,version=as(packageVersion("changepoint.cov"),'character'),date=date(),method=NULL))

#' @describeIn cptCovariance Summarises the cptCovariance object
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("summary","cptCovariance",function(object){
		  cat('Created using changepoint.cov version',object@version,'\n')
		  cat('Method			    : ',object@method,'\n')
		  cat('Multiple changepoint method : ',numCpts(object),'\n')
		  cat('Minimum segment length      : ',object@msl,'\n')
		  cat('Changepoints		    : ',object@cpts,'\n')
})

#' @describeIn cptCovariance Shows the cptCovariance object
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("show","cptCovariance",function(object){
		  cat('Class, cptCovariance  : Covariance Changepoint Object\n')
		  cat('		      : S4 class containing ',length(attributes(object))-1,' slots with names\n')
		  cat('	               ',names(attributes(object))[1:(length(attributes(object))-1)],'\n\n')
		  cat('Created on      : ',object@date,'\n\n')
		  cat('Summary(.)      :\n---------------\n')
		  summary(object)
})

if(!isGeneric("data")){
	if(is.function("data")){
		fun <- data
	}else{
		fun <- function(object){
			standardGeneric("data")
		}
	}
	setGeneric("data",fun)
}

#' @describeIn cptCovariance Retrieves the data slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("data","cptCovariance",function(object){
		  object@data
})

if(!isGeneric("cpts")){
	if(is.function("cpts")){
		fun <- cpts
	}else{
		fun <- function(object){
			standardGeneric("cpts")
		}
	}
	setGeneric("cpts",fun)
}

#' @describeIn cptCovariance Retrieves the cpts slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("cpts","cptCovariance",function(object){
		  object@cpts
})

if(!isGeneric("method")){
	if(is.function("method")){
		fun <- method
	}else{
		fun <- function(object){
			standardGeneric("method")
		}
	}
	setGeneric("method",fun)
}

#' @describeIn cptCovariance Retrieves the method slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("method","cptCovariance",function(object){
		  object@method
})

if(!isGeneric("msl")){
	if(is.function("msl")){
		fun <- msl
	}else{
		fun <- function(object){
			standardGeneric("msl")
		}
	}
	setGeneric("msl",fun)
}

#' @describeIn cptCovariance Retrieves the msl slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("msl","cptCovariance",function(object){
		  object@msl
})

if(!isGeneric("numCpts")){
	if(is.function("numCpts")){
		fun <- numCpts
	}else{
		fun <- function(object){
			standardGeneric("numCpts")
		}
	}
	setGeneric("numCpts",fun)
}

#' @describeIn cptCovariance Retrieves the numCpts slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("numCpts","cptCovariance",function(object){
		  if(is.numeric(object@numCpts)){
			  return(paste0("Manual - ",object@numCpts," changepoints"))
		  }else{
			  return(object@numCpts)
		  }
})

if(!isGeneric("threshold")){
	if(is.function("threshold")){
		fun <- threshold
	}else{
		fun <- function(object){
			standardGeneric("threshold")
		}
	}
	setGeneric("threshold",fun)
}

#' @describeIn cptCovariance Retrieves the threshold slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("threshold","cptCovariance",function(object){
		  object@threshold
})

if(!isGeneric("thresholdValue")){
	if(is.function("thresholdValue")){
		fun <- thresholdValue
	}else{
		fun <- function(object){
			standardGeneric("thresholdValue")
		}
	}
	setGeneric("thresholdValue",fun)
}

#' @describeIn cptCovariance Retrieves the thresholdValue slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("thresholdValue","cptCovariance",function(object){
		  object@thresholdValue
})

if(!isGeneric("cptsSig")){
	if(is.function("cptsSig")){
		fun <- cptsSig
	}else{
		fun <- function(object){
			standardGeneric("cptsSig")
		}
	}
	setGeneric("cptsSig",fun)
}

#' @describeIn cptCovariance Retrieves the cptsSig slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("cptsSig","cptCovariance",function(object){
		  object@cptsSig
})

if(!isGeneric("subspaceDim")){
	if(is.function("subspaceDim")){
		fun <- subspaceDim
	}else{
		fun <- function(object){
			standardGeneric("subspaceDim")
		}
	}
	setGeneric("subspaceDim",fun)
}

#' @describeIn cptCovariance Retrieves the subspaceDim slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("subspaceDim","cptCovariance",function(object){
		  if(toupper(object@method)!='SUBSPACE'){
			  stop("subspaceDim is only a valid slot for method='Subspace'")
		  }else{
			  return(object@subspaceDim)
		  }
})

if(!isGeneric("nperm")){
	if(is.function("nperm")){
		fun <- nperm
	}else{
		fun <- function(object){
			standardGeneric("nperm")
		}
	}
	setGeneric("nperm",fun)
}

#' @describeIn cptCovariance Retrieves the nperm slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("nperm","cptCovariance",function(object){
		  if(!((toupper(object@method)=='SUBSPACE')&&(toupper(object@threshold)=='PERMTEST'))){
			  stop("nperm is only a valid slot when using the permutation test within method='Subspace'")
		  }else{
			  return(object@nperm)
		  }
})

if(!isGeneric("LRCov")){
	if(is.function("LRCov")){
		fun <- LRCov
	}else{
		fun <- function(object){
			standardGeneric("LRCov")
		}
	}
	setGeneric("LRCov",fun)
}

#' @describeIn cptCovariance Retrieves the LRCov slot 
#'
#' @param object An object of S4 class \code{\link{cptCovariance-class}}
setMethod("LRCov","cptCovariance",function(object){
		  if(toupper(object@method)!='CUSUM'){
			  stop("LRCov is only a valid slot for method='CUSUM'")
		  }else{
			  return(object@LRCov)
		  }
})

