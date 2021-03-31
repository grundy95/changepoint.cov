set.seed(1)
dataAMOC <- wishartDataGeneration(n=100,p=3,tau=50)$data
data2Change <- wishartDataGeneration(n=200,p=3,tau=c(60,120))$data
dataHighP <- wishartDataGeneration(n=100,p=20,tau=50)$data

test_that("Default arguments produce no errors",{
		  expect_s4_class(cptCUSUM(X=dataAMOC),"cptCovariance")
})

test_that("cusumTestStat output is correct",{
		  n <- 100
		  msl <- ncol(dataAMOC)
		  ans <- cusumTestStat(dataAMOC,LRCov='Bartlett',msl=msl)
		  expect_type(ans,"double")
		  expect_equal(length(ans),n)
		  expect_equal(sum(is.na(ans)),2*msl-1)
		  expect_error(cusumTestStat(dataAMOC,LRCov=diag(rep(1,5))),"Dimension of manual LRCov is not compatible with data")
})

test_that("Large p datasets output is correct",{
		  n <- 100
		  msl <- 20
		  expect_error(cusumTestStat(dataHighP,LRCov='Bartlett',msl=msl),"Long run covariance estimation not possible, try a different long run covariance estimator or the Ratio method")
		  expect_error(cusumTestStat(dataHighP,LRCov='Empirical',msl=msl),"Long run covariance estimator is not invertible. This is most likely due to the dimension of the data being too large. Please try the Ratio method")
})

test_that("Data is correct format",{
		  dataAMOCdataFrame <- as.data.frame(dataAMOC)
		  set.seed(1)
  		  dataAMOCunivariate <- rnorm(200)
		  dataAMOCcharacter <- matrix(rep('test',200*20),ncol=20)
		  dataAMOCna <- dataAMOC
		  dataAMOCna[1,1] <- NA
		  dataHighDim <- matrix(rnorm(100*75),ncol=75)
		  dataConstant <- matrix(0,ncol=5,nrow=100)

		  expect_s4_class(cptCUSUM(dataAMOC),"cptCovariance")
		  expect_s4_class(cptCUSUM(dataAMOCdataFrame),"cptCovariance")
		  expect_error(cptCUSUM(dataAMOCunivariate),"Data should be a matrix")
		  expect_error(cptCUSUM(as.matrix(dataAMOCunivariate,ncol=1)),"Univariate time series analysis not supported")
		  expect_error(cptCUSUM(dataAMOCcharacter),"Data must be numeric")
		  expect_error(cptCUSUM(dataAMOCna),"Missing value: NA is not allowed in the data")
		  expect_error(cptCUSUM(dataHighDim),"Dimension of data is too high to allow covariance changepoint detection using available methods. Dimension of data needs to be at most n/2")
		  expect_error(cptCUSUM(dataConstant),"Sample covariance of whole data is singular. This is probably due to constant data",fixed=TRUE)
})

test_that("Threshold type is correct",{
		  expect_s4_class(cptCUSUM(dataAMOC,threshold='Asymptotic'),"cptCovariance")
		  expect_s4_class(cptCUSUM(dataAMOC,threshold='Manual'),"cptCovariance")
		  expect_error(cptCUSUM(dataAMOC,threshold='man'),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,threshold=1),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,threshold='Normal'),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,threshold=Normal))
})

test_that("Number of changepoints is correct",{
		  expect_s4_class(cptCUSUM(dataAMOC,numCpts='AMOC'),"cptCovariance")
		  expect_s4_class(cptCUSUM(data2Change,numCpts='BinSeg'),"cptCovariance")
		  expect_s4_class(cptCUSUM(data2Change,numCpts=2),"cptCovariance")

		  expect_error(cptCUSUM(dataAMOC,numCpts='AMC'),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,numCpts=TRUE),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,numCpts=c(4,5)),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,numCpts=amoc))
})

test_that("Threshold value is correct",{
		  expect_s4_class(cptCUSUM(dataAMOC,threshold='Asymptotic',thresholdValue=0.05),"cptCovariance")
		  expect_s4_class(cptCUSUM(dataAMOC,threshold='Manual',thresholdValue=30),"cptCovariance")

		  expect_error(cptCUSUM(dataAMOC,thresholdValue=-1),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
		  expect_error(cptCUSUM(dataAMOC,threshold='Asymptotic',thresholdValue=2),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
		  expect_error(cptCUSUM(dataAMOC,thresholdValue="two"),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
		  expect_error(cptCUSUM(dataAMOC,thresholdValue=c(0.05,0.1)),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'")
})

test_that("Minimum segment length is appropriate",{
		  expect_s4_class(cptCUSUM(dataAMOC,msl=ncol(dataAMOC)+5),"cptCovariance")

		  expect_error(cptCUSUM(dataAMOC,msl=nrow(dataAMOC)-1),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptCUSUM(dataAMOC,msl=-5),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptCUSUM(dataAMOC,msl=ncol(dataAMOC)-2),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptCUSUM(dataAMOC,msl='minimum'),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptCUSUM(dataAMOC,msl=c(ncol(dataAMOC)+1,ncol(dataAMOC)+2)),"Minimum segment length should be a single integer between p and n/2")
})

test_that("LRCov argument is appropriate",{
		  expect_s4_class(cptCUSUM(dataAMOC,LRCov='Bartlett'),"cptCovariance")
		  expect_s4_class(cptCUSUM(dataAMOC,LRCov='Empirical'),"cptCovariance")
		  expect_s4_class(cptCUSUM(dataAMOC,LRCov=diag(rep(1,6))),"cptCovariance")

		  expect_error(cptCUSUM(dataAMOC,LRCov='Normal'),"LRCov not identified: see ?cptCov for valid entries to LRCov",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,LRCov=TRUE),"LRCov not identified: see ?cptCov for valid entries to LRCov",fixed=TRUE)
		  expect_error(cptCUSUM(dataAMOC,LRCov=bartlett))
		  expect_error(cptCUSUM(dataAMOC,LRCov=diag(rep(1,5))),"Dimension of manual LRCov is not compatible with data")
		  expect_error(cptCUSUM(dataHighP,LRCov='Bartlett'),"Long run covariance estimation not possible, try a different long run covariance estimator or the Ratio method")
})

test_that("errorCheck argument is logical",{
		  expect_s4_class(cptCUSUM(dataAMOC,errorCheck=TRUE),"cptCovariance")
		  expect_s4_class(cptCUSUM(dataAMOC,errorCheck=FALSE),"cptCovariance")

		  expect_error(cptCUSUM(dataAMOC,errorCheck='FALSE'))
		  expect_error(cptCUSUM(dataAMOC,errorCheck='TRUE'))
		  expect_error(cptCUSUM(dataAMOC,errorCheck=yes))
})

test_that("Class argument is logical",{
		  expect_s4_class(cptCUSUM(dataAMOC,Class=TRUE),"cptCovariance")
		  expect_type(cptCUSUM(dataAMOC,Class=FALSE),"integer")

		  expect_error(cptCUSUM(dataAMOC,Class='S4'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptCUSUM(dataAMOC,Class='TRUE'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptCUSUM(dataAMOC,Class=yes))
})
##}}}



