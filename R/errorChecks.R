#' Error checking for Ratio Method
#'
#' DEVELOPER USE ONLY. This function checks the user inputs to make sure they are all valid.
#'
#' @inheritParams cptRatio
ratioErrorChecks <- function(X,threshold,numCpts,thresholdValue,msl,Class){
	#data checks
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

	#threshold checks
	if(!is.character(threshold)){
		stop("Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings")
	}
	threshold <- toupper(threshold)
	if((threshold!="ASYMPTOTIC")&&(threshold!="MANUAL")){
		stop("Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings")
	}

	#numCpts checks
	if(!is.character(numCpts)){
		stop("numCpts not identified: see ?cptCov for valid entries to numCpts. NOTE numCpts should be character strings")
	}
	numCpts <- toupper(numCpts)
	if((numCpts!="AMOC")){
		stop("numCpts not identified: see ?cptCov for valid entries to numCpts. NOTE numCpts should be character strings")
	}

	#thresholdValue checks
	if((!is.numeric(thresholdValue))||(length(thresholdValue)!=1)||(thresholdValue<0)){
		stop("thresholdValue must be a single numeric and positive if threshold='Manual'")
	}
	if((threshold=="ASYMPTOTIC")&&(!((thresholdValue==0)||thresholdValue==0.05))){
		warning("thresholdValue is not used when using threshold='Asymptotic' for the Ratio method")
	}

	#msl checks
	if(!is.numeric(msl)){
		stop("Minimum segment length should be a single integer between p and n/2")
	}
	if((length(msl)!=1)||(msl%%1!=0)||(msl<ncol(X))||(msl>(nrow(X)/2))){
		stop("Minimum segment length should be a single integer between p and n/2")
	}

	#Class checks
	if(!is.logical(Class)){
		stop("Class should be logical, TRUE or FALSE")
	}
}

#' Error checking for CUSUM Method
#'
#' DEVELOPER USE ONLY. This function checks the user inputs to make sure they are all valid
#'
#' @inheritParams cptCUSUM
cusumErrorChecks <- function(X,threshold,numCpts,LRCov,thresholdValue,msl,Class){
	#data checks
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

	#threshold checks
	if(!is.character(threshold)){
		stop("Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings")
	}
	threshold <- toupper(threshold)
	if((threshold!="ASYMPTOTIC")&&(threshold!="MANUAL")){
		stop("Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings")
	}

	#numCpts checks
	if(!is.character(numCpts)){
		stop("numCpts not identified: see ?cptCov for valid entries to numCpts. NOTE numCpts should be character strings")
	}
	numCpts <- toupper(numCpts)
	if((numCpts!="AMOC")){
		stop("numCpts not identified: see ?cptCov for valid entries to numCpts. NOTE numCpts should be character strings")
	}

	#LRCov checks
	if(!is.character(LRCov)){
		stop("LRCov not identified: see ?cptCov for valid entries to LRCov. NOTE LRCov should be a character string")
	}
	LRCov <- toupper(LRCov)
	if((LRCov!="BARTLETT")&&(LRCov!="EMPIRICAL")){
		stop("LRCov not identified: see ?cptCov for valid entries to LRCov. NOTE LRCov should be a character string")
	}

	#thresholdValue checks
	if((!is.numeric(thresholdValue))||(length(thresholdValue)!=1)||(thresholdValue<0)){
		stop("thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
	}
	if((threshold=="ASYMPTOTIC")&&(thresholdValue>1)){
		stop("When using threshold='Asymptotic', the thresholdValue is the (1-thresholdValue)-quantile of the standard Normal distribution. Please enter a value between 0 and 1")
	}

	#msl checks
	if(!is.numeric(msl)){
		stop("Minimum segment length should be a single integer between p and n/2")
	}
	if((length(msl)!=1)||(msl%%1!=0)||(msl<ncol(X))||(msl>(nrow(X)/2))){
		stop("Minimum segment length should be a single integer between p and n/2")
	}

	#Class checks
	if(!is.logical(Class)){
		stop("Class should be logical, TRUE or FALSE")
	}
}

#' Error checking for Subspace Method
#'
#' DEVELOPER USE ONLY. This function checks the user inputs to make sure they are all valid
#'
#' @inheritParams cptSubspace
subspaceErrorChecks <- function(X,q,threshold,numCpts,thresholdValue,msl,nperm,Class){
	#data checks
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

	#q checks
	if(!is.numeric(q)){
		stop("Subspace dimension should be a single positive integer")
	}
	if((length(q)!=1)||(q%%1!=0)||(q<1)){
		stop("Subspace dimension should be a single positive integer")
	}
	if(q>ncol(X)){
		stop("Subspace dimension, q, should be smaller than time series dimension, p.")
	}

	#threshold checks
	if(!is.character(threshold)){
		stop("Threshold not identified: see ?cptSubspace for valid entries to threshold. NOTE thresholds should be character strings")
	}
	threshold <- toupper(threshold)
	if((threshold!="PERMTEST")&&(threshold!="MANUAL")){
		stop("Threshold not identified: see ?cptSubspace for valid entries to threshold. NOTE thresholds should be character strings")
	}

	#numCpts checks
	if(!is.character(numCpts)){
		stop("numCpts not identified: see ?cptSubspace for valid entries to numCpts. NOTE numCpts should be character strings")
	}
	numCpts <- toupper(numCpts)
	if((numCpts!="AMOC")){
		stop("numCpts not identified: see ?cptSubspace for valid entries to numCpts. NOTE numCpts should be character strings")
	}

	#thresholdValue checks
	if((!is.numeric(thresholdValue))||(length(thresholdValue)!=1)||(thresholdValue<0)){
		stop("thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='PermTest'")
	}
	if((threshold=="PERMTEST")&&(thresholdValue>1)){
		stop("When using threshold='PermTest', the threshold is the significance of the permutation test. Please enter a value between 0 and 1")
	}

	#msl checks
	if(!is.numeric(msl)){
		stop("Minimum segment length should be a single integer between p and n/2")
	}
	if((length(msl)!=1)||(msl%%1!=0)||(msl<ncol(X))||(msl>(nrow(X)/2))){
		stop("Minimum segment length should be a single integer between p and n/2")
	}

	#nperm checks
	if(!is.numeric(nperm)){
		stop("Number of permutations should be a single positive integer")
	}
	if((threshold!="PERMTEST")&&(nperm!=200)){
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

	#Class checks
	if(!is.logical(Class)){
		stop("Class should be logical, TRUE or FALSE")
	}
}





	





