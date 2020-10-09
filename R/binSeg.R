#' Binary Segmentation Implementation
#' 
#' This function estimates multiple covariance changepoints. NOTE no error checking is performed. This function is for developer use only.
#'
#' @param X Data matrix of dimension n by p.
#' @param method Choice of "Ratio", "CUSUM" or "Subspace".
#' @param msl Minimum segment length
#' @param threshold Threshold choice. Only relevant for using the permutation test within the Subspace method. If method="Subspace" and threshold="PermTest" then the permutation test is used to generate a threshold at each round of the binary segmentation procedure.
#' @param thresholdValue The threshold value that needs determines if a changepoint is significant
#' @param m User specified number of changepoints. Default is m=-1 meaning the thresholdValue is used as a stopping criteria
#' @param LRCov The long-run covariance estimator to be used for CUSUM method. Currently, only "Bartlett" and "Empirical" are supported.
#' @param q Dimension of the latent subspace for Subspace method.
#' @param nperm Only required for Subspace method when threshold="PermTest". Number of permutations to use in the permutation test.
#' @param alpha Only required for Subspace method when threshold="PermTest". Significance level of the permutation test.
#'
#' @return A data frame containing the identified changepoints along with their test statistic and threshold they exceeded.


binSeg <- function(X,method,msl,threshold,thresholdValue,m=-1,LRCov='Bartlett',q=1,nperm=200,alpha=0.05){
	threshold <- toupper(threshold)
	method <- toupper(method)
	n <- nrow(X)
	if(m<0){
		k <- n
	}else{
		k <- m
	}
	cptsSeg <- matrix(c(0,0,0,0),ncol=4)
	cptsSig <- data.frame('cpts'=NaN,'T'=NaN,'thresholdValue'=NaN)
	cpts <- c(0,n)
	while(k>0){
		cpts <- sort(cpts)
		bestSeg <- c(0,n,-1,-1e100)
		for(i in 1:(length(cpts)-1)){
			if(cpts[i+1]-cpts[i]>2*msl){
				if(!(((cpts[i]+1)%in%cptsSeg[,1])&&cpts[i+1]%in%cptsSeg[,2])){
					if(method=='RATIO'){
						testStat <- ratioTestStat(X[(cpts[i]+1):cpts[i+1],],msl)
						cptsSeg <- rbind(cptsSeg,c(cpts[i]+1,cpts[i+1],which.max(testStat)+cpts[i],max(testStat,na.rm=TRUE)))
					}else if(method=='CUSUM'){
						testStat <- cusumTestStat(X[(cpts[i]+1):cpts[i+1],],msl=msl,LRCov=toupper(LRCov))
						cptsSeg <- rbind(cptsSeg,c(cpts[i]+1,cpts[i+1],which.max(testStat)+cpts[i],max(testStat,na.rm=TRUE)))
					}else if(method=='SUBSPACE'){
						 testStat <- subspaceTestStat(X[(cpts[i]+1):cpts[i+1],],msl=msl,q=q)
       						 cptsSeg <- rbind(cptsSeg,c(cpts[i]+1,cpts[i+1],which.max(testStat)+cpts[i],mean(testStat,na.rm=TRUE)))
					}
				}
				if(cptsSeg[cptsSeg[,1]==(cpts[i]+1)&cptsSeg[,2]==cpts[i+1],4]>bestSeg[4]){
					bestSeg <- cptsSeg[cptsSeg[,1]==(cpts[i]+1)&cptsSeg[,2]==cpts[i+1],]
				}
			}
		}
		if(bestSeg[3]<0){break}
		if(m>=0){
			cpts <- c(cpts,bestSeg[3])
			cptsSig <- rbind(cptsSig,c(bestSeg[3:4],NaN))
			k <- k-1
		}else{
		       	if(method=='SUBSPACE'&&threshold=='PERMTEST'){
				thresholdValue <- permutationTest(X[bestSeg[1]:bestSeg[2],],q=q,msl=msl,alpha=alpha,nperm=nperm)
			}
			if(bestSeg[4]<thresholdValue){
				break
			}else{
				cpts <- c(cpts,bestSeg[3])
				cptsSig <- rbind(cptsSig,c(bestSeg[3:4],thresholdValue))
				k <- k-1
			}
		}
	}
	if(!(k==n)){
		cptsSig <- cptsSig[-1,]
	}
	return(cptsSig)
}




