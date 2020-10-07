#' Data Generation - Subspace
#'
#' Allows for easy generation of data that lies in a latent subspace with changepoints in the subspace.
#'
#' @param n Number of time points
#' @param p Dimension of the time series
#' @param q Dimension of the latent subspace
#' @param tau Vector of changepoint locations
#' @param nvar Variance of measurement error (assumed to be Normal)
#' @param svar Variance within the subspace i.e. the variance of points in directions of the basis of the subspace (assumed to be Normal)
#' @param changeSize Either scalar for change size of all the changepoints or a vector of equal length to the number of changepoints indicating the change size for each changepoint. All change sizes must be between 0 and sqrt(q).
#' 
#' @return List with elements:
#' \itemize{
#' 	\item data Matrix of generated data with dimension n by p
#' 	\item cpts Location of the changepoints
#' 	\item changeSize Scalar or vector of change sizes
#' 	\item subspace List of the bases of the subspace for each segment
#' }
#'
#' @importFrom stats rnorm 
subspaceDataGeneration <- function(n,p,q,tau=c(1,n),nvar=0.05,svar=1,changeSize=0.5*sqrt(q)){
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
	if(any(changeSize>sqrt(q))|any(changeSize<0)){
		stop('Change sizes must be between 0 and sqrt(q)')
	}
	if(q>(p/2)){
		qTemp <- p-q
	}else{
		qTemp <- q
	}
	w <- list()
	wTemp <- far::orthonormalization(matrix(rnorm(p*2*qTemp),ncol=2*qTemp),basis=FALSE)
	w[[1]] <- wTemp[,1:qTemp]
	if(m>0){
		w[[2]] <- sqrt(1-min(1,(changeSize[1]^2)/qTemp))*w[[1]]+(changeSize[1]/sqrt(qTemp))*wTemp[,(qTemp+1):(2*qTemp)]
	}
	if(m>1){
		for(i in 3:(m+1)){
			wTemp <- pracma::nullspace(t(w[[i-1]]))[,1:qTemp]
			w[[i]] <- sqrt(1-min(1,(changeSize[i-1]^2)/qTemp))*w[[i-1]]+(changeSize[i-1]/sqrt(qTemp))*wTemp
		}
	}
	if(q>(p/2)){
		for(i in 1:(m+1)){
			w[[i]] <- pracma::nullspace(t(w[[i]]))
		}
	}
	X <- matrix(0,ncol=p,nrow=n)
	for(j in 1:tau[2]){
		if(q==1){
			X[j,] <- w[[1]]*rnorm(q,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
		}else{
			X[j,] <- w[[1]]%*%rnorm(q,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
		}
	}
	if(length(tau)>2){
		for(i in 2:(length(tau)-1)){
			for(j in (tau[i]+1):tau[i+1]){
				if(q==1){
					X[j,] <- w[[i]]*rnorm(q,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
				}else{
					X[j,] <- w[[i]]%*%rnorm(q,0,sqrt(svar))+rnorm(p,0,sqrt(nvar))
				}
			}
		}
	}
	return(list('data'=X,'cpts'=tau,'changeSize'=changeSize,'subspaces'=w))
}

#' Data Generation - Wishart
#'
#' Allows for easy data generation of data containing covariance changepoints as described in Ryan(2020)
#'
#' @param n Number of time points
#' @param p Dimension of the time series
#' @param tau Vector of changepoint locations
#' @param Sigma List of segment covariances. Defaults to automatic generation.
#' @param shape Shape parameter for generation of eigenvalues in transition matrix
#' @param scale Scale parameter for generation of eigenvalues in transition matrix

#' @return List with elements:
#' \itemize{
#' 	\item data Matrix of generated data with dimension n by p
#' 	\item cpts Location of the changepoints
#'	\item Sigma List of covariance matrices of each segment
#' }
#'
#' @importFrom stats rWishart 
wishartDataGeneration <- function(n,p,tau=c(1,n),Sigma=list(NA),shape=5,scale=1/5){
	if(tau[length(tau)]!=n){
		tau <- c(tau,n)
	}
	if(tau[1]!=1){
		tau <- c(1,tau)
	}
	m <- length(tau)-2
	if((length(Sigma)-1!=m)&&(!(is.na(Sigma[[1]])))){stop("Need 1 covariance matrix per segment")}
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
#' @importFrom stats prcomp rWishart rgamma
generateTransMatrix <- function(p,shape=5,scale=1/5){
	A <- prcomp(rWishart(1,p,diag(p))[,,1])[[2]]
	eigs <- rgamma(n=p,shape=shape,scale=scale)
	return(t(A)%*%diag(eigs)%*%A)
}
