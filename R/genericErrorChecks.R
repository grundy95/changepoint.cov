#' data error checking
#'
#' DEVELOPER USE ONLY. This function checks the user input for the data for all methods.
#'
#' @inheritParams cptCov
dataErrorChecks <- function(X){
	if(!is.numeric(X)){
		stop("Data must be numeric")
	}
	if(any(is.na(X))){
		stop("Missing value: NA is not allowed in the data")
	}
	if(!is.matrix(X)){
		stop("Data should be a matrix")
	}
	if(ncol(X)==1){
		stop("Univariate time series analysis not supported")
	}
	if(ncol(X)>(nrow(X)/2)){
		   stop("Dimension of data is too high to allow covariance changepoint detection using available methods. Dimension of data needs to be at most n/2")
	}
}

#' threshold error checking
#'
#' DEVELOPER USE ONLY. This function checks the user input for the threshold and thresholdValue.
#'
#' @param threshold Type of threshold to be used
#' @param thresholdValue Numeric value of threshold or quantile of asymptotic distribution/permutation test
#' @param method function name for relevant help file
thresholdErrorChecks <- function(threshold,thresholdValue,method){
	if(method=='cptCov'){
		if(!is.character(threshold)){
			stop("Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings")
		}
		threshold <- toupper(threshold)
		if((threshold!="ASYMPTOTIC")&&(threshold!="MANUAL")){
			stop("Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings")
		}
		if((!is.numeric(thresholdValue))||(length(thresholdValue)!=1)||(thresholdValue<0)){
			stop("thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
		}
		if((threshold=="ASYMPTOTIC")&&(thresholdValue>1)){
			stop("thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
		}
	}else if(method=='cptSubspace'){
		if(!is.character(threshold)){
			stop("Threshold not identified: see ?cptSubspace for valid entries to threshold. NOTE thresholds should be character strings")
		}
		threshold <- toupper(threshold)
		if((threshold!="PERMTEST")&&(threshold!="MANUAL")){
			stop("Threshold not identified: see ?cptSubspace for valid entries to threshold. NOTE thresholds should be character strings")
		}
		if((!is.numeric(thresholdValue))||(length(thresholdValue)!=1)||(thresholdValue<0)){
			stop("thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='PermTest'")
		}
		if((threshold=="PERMTEST")&&(thresholdValue>1)){
			stop("When using threshold='PermTest', the threshold is the significance of the permutation test. Please enter a value between 0 and 1")
		}
	}
}

#' numCpts error checking
#'
#' DEVELOPER USE ONLY. This function checks the user input for the numCpts argument for all methods.
#'
#' @inheritParams cptCov
numCptsErrorChecks <- function(numCpts,method){
	if(!(is.character(numCpts)||is.numeric(numCpts))){
		stop(paste0("numCpts not identified: see ?",method," for valid entries to numCpts"))
	}
	if(is.numeric(numCpts)){
		if((length(numCpts)!=1)||(numCpts%%1!=0)||(numCpts<1)){
			stop(paste0("numCpts not identified: see ?",method," for valid entries to numCpts"))
		}
	}else{
		numCpts <- toupper(numCpts)
		if(!((numCpts=="AMOC")||(numCpts=='BINSEG'))){
			stop(paste0("numCpts not identified: see ?",method," for valid entries to numCpts"))
		}
	}
}

#' msl error checking
#'
#' DEVELOPER USE ONLY. This function checks the user input for the msl argument for all methods.
#'
#' @inheritParams cptCov
mslErrorChecks <- function(X,msl){
	if(!is.numeric(msl)){
		stop("Minimum segment length should be a single integer between p and n/2")
	}
	if((length(msl)!=1)||(msl%%1!=0)||(msl<ncol(X))||(msl>(nrow(X)/2))){
		stop("Minimum segment length should be a single integer between p and n/2")
	}
}

#' Class error checking
#'
#' DEVELOPER USE ONLY. This function checks the user input for the class argument for all methods.
#'
#' @inheritParams cptCov
classErrorChecks <- function(Class){
	if(!is.logical(Class)){
		stop("Class should be logical, TRUE or FALSE")
	}
}
