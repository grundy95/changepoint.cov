#' Test statistic for Ratio method
#'
#' Calculates the Ratio test statistic for all potential changepoint locations within the time series. See \code{\link{cptRatio}} for more details.
#'
#' See \code{\link{cptRatio}}.
#'
#' @param X Data matrix of dimension n by p.
#' @param msl A numeric giving the minimum segment length between changepoints. NOTE this should be greater than or equal to p.
#'
#' @return A numeric vector containing the test statistic at each potential changepoint location.
#'
#' @references
#' \insertRef{Ryan2020}{changepoint.cov}
#'
#' @seealso \code{\link{cptRatio}}, \code{\link{cptCov}}
#'
#' @examples
#' set.seed(1)
#' data <- wishartDataGeneration(n=100,p=3,tau=50)$data
#' ans <- ratioTestStat(X=data,msl=20)
#' which.max(ans)
#'
#' @importFrom rlang !!!
#' @export
ratioTestStat <- function(X,msl){
	n <- nrow(X)
	p <- ncol(X)
	X <- purrr::map(1:n,~X[.,])
	calculateRatioDistance <- ratioDistanceCalculator(X)	
	T <- purrr::map(msl:(n-msl),calculateRatioDistance)
	gammas <- purrr::map(msl:(n-msl),~c(p/.,p/(n-.)))
	trace <- purrr::map(gammas,~rlang::exec(calculateExpectedTrace,!!!.))
	bias <- purrr::map(gammas,~rlang::exec(asymptoticBias,!!!.))
	variance <- purrr::map(gammas,~rlang::exec(asymptoticVariance,!!!.))

	traceErrors <- purrr::map_lgl(trace,is.null)
	biasErrors <- purrr::map_lgl(bias,is.na)
	varianceErrors <- purrr::map_lgl(variance,is.na)
	testErrors <- purrr::pmap_lgl(list(traceErrors,biasErrors,varianceErrors),~(..1||..2||..3))
	if(any(testErrors==TRUE)){
		warning('msl is too short making the test statistic incalculable for certain boundary changepoint locations. Method has continued with the exclusion of the incalculable potential changepoint locations.')
	}	

	trace <- purrr::map_if(trace,testErrors,~NA_real_)
	trace <- purrr::map(trace,~.[[1]])
	bias <- purrr::map_if(bias,testErrors,~NA_real_)
	variance <- purrr::map_if(variance,testErrors,~NA_real_)
	testStat <- purrr::pmap_dbl(list(T,trace,bias,variance),~(..4^(-0.25))*(..1-p*..2-..3))
	return(c(rep(NA_real_,(msl-1)),testStat,rep(NA_real_,msl)))
}

#' Ratio distance calculator
#'
#' Creates a function that takes a potential changepoint location and returns the un-normalized test statistic defined in \insertCite{Ryan2020;textual}{changepoint.cov}.
#'
#' @param X List of data where each slot is a time point
#'
#' @return A function used to calculate un-normalized test statistic
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
ratioDistanceCalculator <- function(X){
	A <- purrr::map(X,~.%*%t(.))
	A <- purrr::accumulate(A,`+`)
	n <- length(X)
	function(tau){
		Sigma1 <- (1/tau)*A[[tau]]
		Sigma2 <- (1/(n-tau))*(A[[n]]-A[[tau]])
		eig <- geigen::geigen(Sigma2,Sigma1,symmetric=TRUE,only.values=TRUE)$values
		return(rlang::exec(sum,(1-eig)^2)+rlang::exec(sum,(1-(1/eig))^2))
	}
}

#' Expected Trace Calulator
#' 
#' Calculates the expected trace used in the Ratio test statistic
#'
#' @param gamma1 p/n1 where n1 is the number of time points before potential changepoint
#'
#' @param gamma2 p/n2 where n2 is the number of time points after potential changepoint
#'
#' @return A list where the first element is the integral and second element is the error.
#'
#' @importFrom stats integrate
#' @importFrom rlang !!!
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
calculateExpectedTrace <- function(gamma1,gamma2){
	asymptoticPdf <- fisherESD(gamma1,gamma2)
	integrand <- functionProduct(function(x){(1-x)^2+(1-1/x)^2},asymptoticPdf)
	asymptoticSupport <- fisherSupport(gamma1,gamma2)
	safeIntegral <- purrr::safely(integrate)
	integral <- rlang::exec(safeIntegral,integrand,`!!!`(asymptoticSupport))[[1]]
	return(integral)
}

#' Empirical spectral distribution of a Fisher matrix
#'
#' Creates a function which returns the ESD of a Fisher matrix with set gammas for some x within the support of the ESD. Used within Ratio method.
#'
#' @inheritParams calculateExpectedTrace
#'
#' @return A function that calculates the ESD of a Fisher matrix for given gammas
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
fisherESD <- function(gamma1,gamma2){
	h <- sqrt(gamma1+gamma2-gamma1*gamma2)
	a <- ((1-h)^2)/((1-gamma2)^2)
	b <- ((1+h)^2)/((1-gamma2)^2)
	safeSqrt <- purrr::quietly(sqrt)
	function(x){
		result <- rep(0,length(x))
		result[x>a&x<b] <- ((1-gamma2)/(2*pi*x*(gamma1+gamma2*x)))*safeSqrt((b-x)*(x-a))[[1]]
		return(result)
	}
}

#' Calculates Fisher support
#'
#' Returns the support for calculating the Fisher ESD within Ratio method.
#'
#' @inheritParams calculateExpectedTrace
#'
#' @return Support for Fisher ESD
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
fisherSupport <- function(gamma1,gamma2){
	h <- sqrt(gamma1+gamma2-gamma1*gamma2)
	a <- ((1-h)^2)/((1-gamma2)^2)
	b <- ((1+h)^2)/((1-gamma2)^2)
	return(c(a,b))
}

#' Product of two functions
#'
#' Multiplies outputs of two functions 
#'
#' @param f1 Function 1
#' @param f2 Function 2
#'
#' @return Function which is product of functions
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
functionProduct <- function(f1,f2){
	return(function(x){return(f1(x)*f2(x))})
}

#' Asymptotic Bias of Ratio Test Staistic
#'
#' Calculates the asymptotic bias of the Ratio test statistic used to generate the normalized test statistic
#'
#' @inheritParams calculateExpectedTrace
#'
#' @return A numeric of the asymptotic bias
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
asymptoticBias <- function(gamma1,gamma2){
	h <- sqrt(gamma1+gamma2-gamma1*gamma2)
	K31 <- (h^2)/((1-gamma2)^4)
	K21 <- (2*h*(1+h^2))/((1-gamma2)^4)-(2*h)/((1-gamma2)^2)
	K32 <- (h^2)/((1-gamma1)^4)
	K22 <- (2*h*(1+h^2))/((1-gamma1)^4)-(2*h)/((1-gamma1)^2)
	return(2*K31*(1-(gamma2^2)/(h^2))+(2*K21*gamma2)/(h)+2*K32*(1-(gamma1^2)/(h^2))+(2*K22*gamma1)/h)
}

#' Asymptotic Variance of Ratio Test Statistic
#'
#' Calculates the asymptotic variance of the Ratio test statistic used to generate the normalized test statistic
#'
#' @inheritParams calculateExpectedTrace
#'
#' @return A numeric of the asymptotic variance
#'
#' @seealso \code{\link{ratioTestStat}}
#'
#' @keywords internal
asymptoticVariance <- function(gamma1,gamma2){
	h <- sqrt(gamma1+gamma2-gamma1*gamma2)
	K31 <- (h^2)/((1-gamma2)^4)
	K21 <- (2*h*(1+h^2))/((1-gamma2)^4)-(2*h)/((1-gamma2)^2)
	K32 <- (h^2)/((1-gamma1)^4)
	K22 <- (2*h*(1+h^2))/((1-gamma1)^4)-(2*h)/((1-gamma1)^2)
	J1 <- -2*(1-gamma2)^2
	J2 <- (1-gamma2)^4
	return(2*(K21^2+2*(K31^2)+K22^2+2*(K32^2)+(J1*K21)/(h)+(J1*K21)/(h*(h^2-1))-(J1*K31*(h^2+1))/(h^2)-(J1*K31)/(h^2*(h^2-1))+(J2*K21*2*h)/((h^2-1)^3)+(J2*K31)/(h^2)+(J2*K31*(1-3*h^2))/((h^2)*(h^2-1)^3)))
}



