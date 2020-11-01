#' Binary Segmentation Implementation
#' 
#' This function estimates multiple covariance changepoints. NOTE no argument error checking is performed. This function is for developer use only.
#'
#' The function works by storing a matrix with columns corresponding to the start point of segment; end point of segment; most likely changepoint position within segment; and test statistic for the segment. This way as more and more segments are detected then we don't repeat the calculation of the test statistic. The while loop allows for a set number of changepoints to be calculated or if this is unknown then it will keep looping until no more changepoints are possible or the maximum test statistic across the segments doesn't exceed the threshold. For the Permutation Test the Binary Segmentation only performs a Permutation Test on the segment with the maximum test statistic.
#'
#' @param X Data matrix of dimension n by p.
#' @param method Choice of "Ratio", "CUSUM" or "Subspace".
#' @param msl Minimum segment length. NOTE must be greater than p/2.
#' @param threshold Threshold choice. Only relevant for using the permutation test within the Subspace method. If method="Subspace" and threshold="PermTest" then the permutation test is used to generate a threshold at each round of the binary segmentation procedure.
#' @param thresholdValue The threshold value that needs determines if a changepoint is significant.
#' @param m User specified number of changepoints. Default is m=-1 meaning the thresholdValue is used as a stopping criteria.
#' @param LRCov The long-run covariance estimator to be used for CUSUM method. Currently, only "Bartlett"; "Empirical" or manual long-run covariance estimate are supported. For a manual estimate the matrix containing the estimate should be inputted - no error checks are performed.
#' @param subspaceDim Dimension of the latent subspace for Subspace method.
#' @param nperm Only required for Subspace method when threshold="PermTest". Number of permutations to use in the permutation test.
#' @param alpha Only required for Subspace method when threshold="PermTest". Significance level of the permutation test.
#'
#' @return A data frame containing the identified changepoints along with their test statistic and threshold they exceeded.
#'
#' @keywords internal


binSeg <- function(X,method,msl,threshold='Manual',thresholdValue=0,m=-1,LRCov='Bartlett',subspaceDim=1,nperm=200,alpha=0.05){
	threshold <- toupper(threshold)
	method <- toupper(method)
	if(is.character(LRCov)){
		LRCov <- toupper(LRCov)
	}
	n <- nrow(X)
	if(msl>(n/2)){
		stop("Minimum segment length should be a single integer between p and n/2")
	}
	if(m<0){
		k <- n
	}else{
		k <- m
	}
	cptsSeg <- matrix(c(0,0,0,0),ncol=4)
	cptsSig <- data.frame('cpts'=NA_integer_,'T'=NA_real_,'thresholdValue'=NA_real_,'significant'=FALSE)
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
						testStat <- cusumTestStat(X[(cpts[i]+1):cpts[i+1],],msl=msl,LRCov=LRCov)
						cptsSeg <- rbind(cptsSeg,c(cpts[i]+1,cpts[i+1],which.max(testStat)+cpts[i],max(testStat,na.rm=TRUE)))
					}else if(method=='SUBSPACE'){
						 testStat <- subspaceTestStat(X[(cpts[i]+1):cpts[i+1],],msl=msl,subspaceDim=subspaceDim)
       						 cptsSeg <- rbind(cptsSeg,c(cpts[i]+1,cpts[i+1],which.max(testStat)+cpts[i],mean(testStat,na.rm=TRUE)))
					}
				}
				if(cptsSeg[cptsSeg[,1]==(cpts[i]+1)&cptsSeg[,2]==cpts[i+1],4]>bestSeg[4]){
					bestSeg <- cptsSeg[cptsSeg[,1]==(cpts[i]+1)&cptsSeg[,2]==cpts[i+1],]
				}
			}
		}
		if(bestSeg[3]<0){
			if(m>=0){
				warning(paste0("Cannot allocate ",m," changepoints due to minimum segment length restrictions"))
			}
			break
		}
		if(m>=0){
			cpts <- c(cpts,bestSeg[3])
			cptsSig <- rbind(cptsSig,data.frame('cpts'=bestSeg[3],'T'=bestSeg[4],'thresholdValue'=NA_real_,'significant'=TRUE))
			k <- k-1
		}else{
		       	if(method=='SUBSPACE'&&threshold=='PERMTEST'){
				thresholdValue <- permutationTest(X[bestSeg[1]:bestSeg[2],],subspaceDim=subspaceDim,msl=msl,alpha=alpha,nperm=nperm)
			}
			if(bestSeg[4]<thresholdValue){
				cptsSig <- rbind(cptsSig,data.frame('cpts'=bestSeg[3],'T'=bestSeg[4],'thresholdValue'=thresholdValue,'significant'=FALSE))
				break
			}else{
				cpts <- c(cpts,bestSeg[3])
				cptsSig <- rbind(cptsSig,data.frame('cpts'=bestSeg[3],'T'=bestSeg[4],'thresholdValue'=thresholdValue,'significant'=TRUE))
				k <- k-1
			}
		}
	}
	return(cptsSig[-1,])
}




