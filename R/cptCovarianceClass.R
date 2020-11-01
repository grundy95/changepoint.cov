#' An S4 class for a covariance changepoint object
#'
#' List of the basic methods, retrieval functions and slots
#' 
#' @slot data An n by p matrix of the data.
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
#' out <- new('cptCovariance',data=matrix(rnorm(300),ncol=3),
#'			cpts=c(50,100),
#'			method='Ratio',
#'			numCpts='AMOC',
#'			cptsSig=data.frame('cpts'=50,'T'=33.3,'thresholdValue'=30,significant=TRUE),
#'			threshold='Manual',
#'			thresholdValue=30,
#'			msl=20)
#' summary(out)
#' show(out)
#' plot(out)
#' cptsSig(out)
#' covEst(out)
#'
#' @import methods
#' @export
setClass("cptCovariance",slots=list(data='matrix',cpts='numeric',method='character',msl='numeric',numCpts='ANY',threshold='character',thresholdValue='numeric',cptsSig='data.frame',subspaceDim='numeric',nperm='numeric',LRCov='ANY',covEst='list',subspaceEst='list',date='character',version='character'),prototype=list(subspaceDim=NA_real_,nperm=NA_real_,LRCov=NA_character_,covEst=list(NA_real_),subspaceEst=list(NA_real_),version=as(packageVersion("changepoint.cov"),'character'),date=date(),method=NULL))

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

if(!isGeneric("cptsSig")){
	if(is.function("cptsSig")){
		fun <- cptsSig
	}else{
		fun <- function(x){
			standardGeneric("cptsSig")
		}
	}
	setGeneric("cptsSig",fun)
}

#' @describeIn cptCovariance Returns a data frame containing the changepoints, the associated test statistic and threshold.
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases cptsSig
#'
#' @export
setMethod("cptsSig","cptCovariance",function(x){
		  return(x@cptsSig)
})


#' @describeIn cptCovariance Plotting method for cptCovariance object. Returns a \code{\link[ggplot2]{ggplot}} object which can be manipulated as required
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod("plot","cptCovariance",function(x){
		  p <- ncol(x@data)
		  covs <- covEst(x)

		  segCovs <- data.frame('Segment'=1,'p1'=rep(1:p,each=p),'p2'=rep(1:p,p),'Value'=as.vector(covs[[1]]))
		  if(length(cpts(x))>1){
			  for(i in 2:length(cpts(x))){
				  segCovs <- rbind(segCovs,data.frame('Segment'=i,'p1'=rep(1:p,each=p),'p2'=rep(1:p,p),'Value'=as.vector(covs[[i]])))
			  }
		  }
		  covPlot <- ggplot(segCovs,aes(x=.data$p1,y=.data$p2,fill=.data$Value))+
			  geom_tile()+
			  facet_grid(.~Segment,labeller='label_both')+
			  scale_fill_viridis()+
			  xlab("p")+ylab("p")
		  return(covPlot)
})




if(!isGeneric("covEst")){
	if(is.function("covEst")){
		fun <- covEst
	}else{
		fun <- function(x){
			standardGeneric("covEst")
		}
	}
	setGeneric("covEst",fun)
}

#' @describeIn cptCovariance Returns covariance estimates for each segment
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases covEst
#'
#' @export
setMethod("covEst","cptCovariance",function(x){
		  X <- x@data
		  covs <- list()
		  cpts <- c(0,cpts(x))
		  for(i in 1:(length(cpts)-1)){
			  covs[[i]] <- cov(X[(cpts[i]+1):cpts[i+1],])
		  }
		  return(covs)
})

if(!isGeneric("subspaceEst")){
	if(is.function("subspaceEst")){
		fun <- subspaceEst
	}else{
		fun <- function(x){
			standardGeneric("subspaceEst")
		}
	}
	setGeneric("subspaceEst",fun)
}

#' @describeIn cptCovariance Returns a basis of the subspace estimates for each segment
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases subspaceEst
#'
#' @export
setMethod("subspaceEst","cptCovariance",function(x){
		  if(method(x)!='Subspace'){
			  stop("Subspace estimation only possible for method='Subspace'")
		  }else{
			  X <- x@data
			  q <- subspaceDim(x)
			  subspace <- list()
			  cpts <- c(0,cpts(x))
			  for(i in 1:(length(cpts)-1)){
				 subspace[[i]] <- eigen(cov(X[(cpts[i]+1):cpts[i+1],]),symmetric=TRUE)$vectors[,1:q]
			  }
		  }
		  return(subspace)
})

if(!isGeneric("data")){
	if(is.function("data")){
		fun <- data
	}else{
		fun <- function(x){
			standardGeneric("data")
		}
	}
	setGeneric("data",fun)
}

#' @describeIn cptCovariance Retrieves the data slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases data
#'
#' @export
setMethod("data","cptCovariance",function(x){
		  x@data
})

if(!isGeneric("cpts")){
	if(is.function("cpts")){
		fun <- cpts
	}else{
		fun <- function(x){
			standardGeneric("cpts")
		}
	}
	setGeneric("cpts",fun)
}

#' @describeIn cptCovariance Retrieves the cpts slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases cpts
#'
#' @export
setMethod("cpts","cptCovariance",function(x){
		  x@cpts
})

if(!isGeneric("method")){
	if(is.function("method")){
		fun <- method
	}else{
		fun <- function(x){
			standardGeneric("method")
		}
	}
	setGeneric("method",fun)
}

#' @describeIn cptCovariance Retrieves the method slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases method
#'
#' @export
setMethod("method","cptCovariance",function(x){
		  x@method
})

if(!isGeneric("msl")){
	if(is.function("msl")){
		fun <- msl
	}else{
		fun <- function(x){
			standardGeneric("msl")
		}
	}
	setGeneric("msl",fun)
}

#' @describeIn cptCovariance Retrieves the msl slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases msl
#'
#' @export
setMethod("msl","cptCovariance",function(x){
		  x@msl
})

if(!isGeneric("numCpts")){
	if(is.function("numCpts")){
		fun <- numCpts
	}else{
		fun <- function(x){
			standardGeneric("numCpts")
		}
	}
	setGeneric("numCpts",fun)
}

#' @describeIn cptCovariance Retrieves the numCpts slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases numCpts
#'
#' @export
setMethod("numCpts","cptCovariance",function(x){
		  if(is.numeric(x@numCpts)){
			  return(paste0("Manual - ",x@numCpts," changepoints"))
		  }else{
			  return(x@numCpts)
		  }
})

if(!isGeneric("threshold")){
	if(is.function("threshold")){
		fun <- threshold
	}else{
		fun <- function(x){
			standardGeneric("threshold")
		}
	}
	setGeneric("threshold",fun)
}

#' @describeIn cptCovariance Retrieves the threshold slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases threshold
#'
#' @export
setMethod("threshold","cptCovariance",function(x){
		  x@threshold
})

if(!isGeneric("thresholdValue")){
	if(is.function("thresholdValue")){
		fun <- thresholdValue
	}else{
		fun <- function(x){
			standardGeneric("thresholdValue")
		}
	}
	setGeneric("thresholdValue",fun)
}

#' @describeIn cptCovariance Retrieves the thresholdValue slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases thresholdValue
#'
#' @export
setMethod("thresholdValue","cptCovariance",function(x){
		  x@thresholdValue
})

if(!isGeneric("subspaceDim")){
	if(is.function("subspaceDim")){
		fun <- subspaceDim
	}else{
		fun <- function(x){
			standardGeneric("subspaceDim")
		}
	}
	setGeneric("subspaceDim",fun)
}

#' @describeIn cptCovariance Retrieves the subspaceDim slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases subspaceDim
#'
#' @export
setMethod("subspaceDim","cptCovariance",function(x){
		  if(toupper(x@method)!='SUBSPACE'){
			  stop("subspaceDim is only a valid slot for method='Subspace'")
		  }else{
			  return(x@subspaceDim)
		  }
})

if(!isGeneric("nperm")){
	if(is.function("nperm")){
		fun <- nperm
	}else{
		fun <- function(x){
			standardGeneric("nperm")
		}
	}
	setGeneric("nperm",fun)
}

#' @describeIn cptCovariance Retrieves the nperm slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases nperm
#'
#' @export
setMethod("nperm","cptCovariance",function(x){
		  if(!((toupper(x@method)=='SUBSPACE')&&(toupper(x@threshold)=='PERMTEST'))){
			  stop("nperm is only a valid slot when using the permutation test within method='Subspace'")
		  }else{
			  return(x@nperm)
		  }
})

if(!isGeneric("LRCov")){
	if(is.function("LRCov")){
		fun <- LRCov
	}else{
		fun <- function(x){
			standardGeneric("LRCov")
		}
	}
	setGeneric("LRCov",fun)
}

#' @describeIn cptCovariance Retrieves the LRCov slot 
#'
#' @param x An object of S4 class \code{\linkS4class{cptCovariance}}
#'
#' @aliases LRCov
#'
#' @export
setMethod("LRCov","cptCovariance",function(x){
		  if(toupper(x@method)!='CUSUM'){
			  stop("LRCov is only a valid slot for method='CUSUM'")
		  }else{
			  return(x@LRCov)
		  }
})

