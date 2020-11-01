#' Data Generation - Subspace
#'
#' Allows for easy generation of data that lies in a latent subspace with changepoints in the subspace.
#'
#' This function allows for the creation of a data matrix that contains changes in subspaces as described in \insertCite{Grundy2020;textual}{changepoint.cov}. The creation of the data set starts by generating a random initial basis for the first segment. Subsequent subspace basis are then created that are the required change size from the previous subspace basis. \code{svar} is then the variance of the points in the directions of the subspace bases and \code{nvar} is the variance of the points in the directions orthogonal to these subspace bases. Hence, \code{nvar} must be less then \code{svar} for the subspace to be visible. The change size is defined in terms of the F-distance between subspace bases as described in \insertCite{Grundy2020;textual}{changepoint.cov}.
#'
#' @param n Number of time points.
#' @param p Dimension of the time series.
#' @param subspaceDim Dimension of the latent subspace.
#' @param tau Vector of changepoint locations.
#' @param nvar Variance of measurement error i.e the variance of the points in the directions of the orthogonal compliment of the subspace (assumed to be Normal).
#' @param svar Variance within the subspace i.e. the variance of points in directions of the basis of the subspace (assumed to be Normal).
#' @param changeSize Either scalar for change size of all the changepoints or a vector (of equal length to the number of changepoints) indicating the change size for each changepoint. All change sizes must be between 0 and sqrt(subspaceDim).
#' 
#' @return List with elements:
#' \itemize{
#' 	\item data Matrix of generated data with dimension n by p.
#' 	\item cpts Location of the changepoints.
#' 	\item changeSize Scalar or vector of change sizes.
#' 	\item subspace List of the bases of the subspace for each segment.
#' }
#'
#' @references 
#' \insertRef{Grundy2020}{changepoint.cov}
#'
#' @seealso \code{\link{cptSubspace}}, \code{\link{wishartDataGeneration}}
#'
#'
#' @examples
#' set.seed(1)
#' dataObject <- subspaceDataGeneration(n=100,p=4,subspaceDim=2,tau=40)
#' dim(dataObject$data)
#' dataObject$cpts
#' dataObject$changeSize
#' dataObject$subspace
#'
#' @importFrom stats rnorm 
#' @export
subspaceDataGeneration <- function(n,p,subspaceDim,tau=c(1,n),nvar=0.05,svar=1,changeSize=0.5*sqrt(subspaceDim)){
	if(tau[1]!=1){
		tau <- c(1,tau)
	}
	if(tau[length(tau)]!=n){
		tau <- c(tau,n)
	}
	m <- length(tau)-2
	if(length(changeSize)!=m&length(changeSize)==1){
		changeSize <- rep(changeSize,m)
	}else if(length(changeSize)!=m&length(changeSize)>1){
		stop('Number of changepoints does not match number of change sizes')
	}
	if(any(changeSize>sqrt(subspaceDim))|any(changeSize<0)){
		stop('Change sizes must be between 0 and sqrt(subspaceDim)')
	}
	if(subspaceDim>(p/2)){
		subspaceDimTemp <- p-subspaceDim
	}else{
		subspaceDimTemp <- subspaceDim
	}
	w <- list()
	wTemp <- far::orthonormalization(matrix(rnorm(p*2*subspaceDimTemp),ncol=2*subspaceDimTemp),basis=FALSE)
	w[[1]] <- wTemp[,1:subspaceDimTemp]
	if(m>0){
		w[[2]] <- sqrt(1-min(1,(changeSize[1]^2)/subspaceDimTemp))*w[[1]]+(changeSize[1]/sqrt(subspaceDimTemp))*wTemp[,(subspaceDimTemp+1):(2*subspaceDimTemp)]
	}
	if(m>1){
		for(i in 3:(m+1)){
			wTemp <- pracma::nullspace(t(w[[i-1]]))[,1:subspaceDimTemp]
			w[[i]] <- sqrt(1-min(1,(changeSize[i-1]^2)/subspaceDimTemp))*w[[i-1]]+(changeSize[i-1]/sqrt(subspaceDimTemp))*wTemp
		}
	}
	if(subspaceDim>(p/2)){
		for(i in 1:(m+1)){
			w[[i]] <- pracma::nullspace(t(w[[i]]))
		}
	}
	X <- matrix(0,ncol=p,nrow=n)
	for(j in 1:tau[2]){
		if(subspaceDim==1){
			X[j,] <- w[[1]]*rnorm(subspaceDim,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
		}else{
			X[j,] <- w[[1]]%*%rnorm(subspaceDim,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
		}
	}
	if(length(tau)>2){
		for(i in 2:(length(tau)-1)){
			for(j in (tau[i]+1):tau[i+1]){
				if(subspaceDim==1){
					X[j,] <- w[[i]]*rnorm(subspaceDim,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
				}else{
					X[j,] <- w[[i]]%*%rnorm(subspaceDim,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
				}
			}
		}
	}
	return(list('data'=X,'cpts'=tau,'changeSize'=changeSize,'subspaces'=w))
}

#' Data Generation - Wishart
#'
#' Allows for easy generation of data containing covariance changepoints as described in \insertCite{Ryan2020;textual}{changepoint.cov}
#'
#' This function generates data that contain changes in covariance. If the covariance for each segment, Sigma, is given then data is generated from a $N(0,Sigma)$ for each segment. If the segment covariance aren't given then they are generated via simulating from Wishart distribution.
#'
#' @param n Number of time points.
#' @param p Dimension of the time series.
#' @param tau Vector of changepoint locations.
#' @param Sigma List of segment covariances. Defaults to automatic generation.
#' @param shape Shape parameter for generation of eigenvalues in transition matrix. Only used with automatic generation of segment covariances.
#' @param scale Scale parameter for generation of eigenvalues in transition matrix. Only used with automatic generation of segment covariances. 

#' @return List with elements:
#' \itemize{
#' 	\item data Matrix of generated data with dimension n by p.
#' 	\item cpts Location of the changepoints.
#'	\item Sigma List of covariance matrices of each segment.
#' }
#'
#' @references
#' \insertRef{Ryan2020}{changepoint.cov}
#'
#' @seealso \code{\link{cptCov}}, \code{\link{subspaceDataGeneration}}
#'
#' @examples
#' set.seed(1)
#' dataObject <- wishartDataGeneration(n=100,p=4,tau=40)
#' dim(dataObject$data)
#' dataObject$cpts
#' dataObject$Sigma
#'
#' @importFrom stats rWishart 
#' @export
wishartDataGeneration <- function(n,p,tau=c(1,n),Sigma=list(NA),shape=5,scale=1/5){
	if(tau[length(tau)]!=n){
		tau <- c(tau,n)
	}
	if(tau[1]!=1){
		tau <- c(1,tau)
	}
	m <- length(tau)-2
	if(((length(Sigma)-1)!=m)&&((any(!is.na(Sigma[[1]]))))){stop("Need 1 covariance matrix per segment")}
	if(is.na(Sigma[[1]])){
		Sigma[[1]] <- (1/p)*rWishart(1,p,diag(p))[,,1]
		if(m>0){
			transMatrix <- append(Sigma,purrr::map(1:m,~generateTransMatrix(p,shape,scale)))
			Sigma <- purrr::accumulate(transMatrix,function(X,Y){B=pracma::sqrtm(Y)$B;return(B%*%Y%*%B)})
		}
	}
	if(m>0){
		cptGaps <- diff(c(0,tau[-1]))
	}else{
		cptGaps <- n
	}
	X <- purrr::map2(cptGaps,Sigma,~MASS::mvrnorm(.x,rep(0,p),.y))
	X <- purrr::reduce(X,rbind)
	return(list(data=X,cpts=tau,Sigma=Sigma))
}
			
#' Transition matrix generation
#'
#' Generates transition matrices for \code{\link{wishartDataGeneration}}
#'
#' @inheritParams wishartDataGeneration
#'
#' @return A p by p transition matrix
#'
#' @seealso \code{\link{wishartDataGeneration}}
#'
#' @importFrom stats prcomp rWishart rgamma
#'
#' @keywords internal
generateTransMatrix <- function(p,shape=5,scale=1/5){
	A <- prcomp(rWishart(1,p,diag(p))[,,1])[[2]]
	eigs <- rgamma(n=p,shape=shape,scale=scale)
	return(t(A)%*%diag(eigs)%*%A)
}
