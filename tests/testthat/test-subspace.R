context('cptSubspace Tests')
library(changepoint.cov)

##{{{ Data Creation
set.seed(1)
dataAMOC <- subspaceDataGeneration(n=100,p=20,subspaceDim=5,tau=50,changeSize=0.5*sqrt(5),nvar=0.05,svar=1)$data
dataNull <- subspaceDataGeneration(n=100,p=20,subspaceDim=5,tau=100,changeSize=0)$data
data2Change <- subspaceDataGeneration(n=100,p=10,subspaceDim=3,tau=c(30,60))$data

##}}}

##{{{ Basic Functionality 
test_that("Default arguments produce no errors",{
		  expect_is(cptSubspace(X=dataAMOC,subspaceDim=5),"cptCovariance")
		  expect_is(cptSubspace(X=dataNull,subspaceDim=5),"cptCovariance")
})

test_that("Manual threshold works for rest of tests",{
		  expect_is(cptSubspace(X=dataAMOC,subspaceDim=5,threshold='Manual',thresholdValue=10),"cptCovariance")
})

test_that("subspaceTestStat output is correct",{
		  n <- 100
		  msl <- ncol(dataAMOC)
		  ans <- subspaceTestStat(dataAMOC,subspaceDim=5,msl=msl)
		  expect_is(ans,'numeric')
		  expect_equal(length(ans),n)
		  expect_equal(sum(is.na(ans)),2*msl-1)
})

test_that("permutationTest output is correct",{
		  n <- 100
		  msl <- ncol(dataAMOC)
		  thresholdValue <- 0.05
		  set.seed(1)
		  ans <- permutationTest(dataAMOC,subspaceDim=5,msl=msl,alpha=0.05,nperm=200)
		  expect_is(ans,'numeric')
		  expect_equal(length(ans),1)
		  expect_true(ans>0)
})
##}}}

##{{{ Error catching tests
test_that("Data is correct format",{
		  dataAMOCdataFrame <- as.data.frame(dataAMOC)
		  set.seed(1)
  		  dataAMOCunivariate <- rnorm(200)
		  dataAMOCcharacter <- matrix(rep('test',200*20),ncol=20)
		  dataAMOCna <- dataAMOC
		  dataAMOCna[1,1] <- NA
		  dataHighDim <- matrix(rnorm(100*75),ncol=75)
		  dataConstant <- matrix(0,ncol=5,nrow=100)
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='Manual',thresholdValue=10),"cptCovariance")
		  expect_is(cptSubspace(dataAMOCdataFrame,subspaceDim=5,threshold='Manual',thresholdValue=10),"cptCovariance")

		  expect_error(cptSubspace(dataAMOCunivariate,subspaceDim=5),"Data should be a matrix")
		  expect_error(cptSubspace(as.matrix(dataAMOCunivariate,ncol=1),subspaceDim=5),"Univariate time series analysis not supported")
		  expect_error(cptSubspace(dataAMOCcharacter,subspaceDim=5),"Data must be numeric")
		  expect_error(cptSubspace(dataAMOCna,subspaceDim=5),"Missing value: NA is not allowed in the data")
		  expect_error(cptSubspace(dataHighDim,subspaceDim=10),"Dimension of data is too high to allow covariance changepoint detection using available methods. Dimension of data needs to be at most n/2")
		  expect_error(cptCUSUM(dataConstant),"Sample covariance of whole data is singular. This is probably due to constant data",fixed=TRUE)
})

test_that("Subspace dimension is correct format",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='Manual',thresholdValue=10),"cptCovariance")

		  expect_error(cptSubspace(dataAMOC,subspaceDim=40),"Subspace dimension, subspaceDim, should be smaller than time series dimension, p.")
		  expect_error(cptSubspace(dataAMOC,subspaceDim="one"),"Subspace dimension should be a single positive integer")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=c(2,3)),"Subspace dimension should be a single positive integer")
})

test_that("Threshold type is correct",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest'),"cptCovariance")
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='Manual'),"cptCovariance")

		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold='man'),"Threshold not identified: see ?cptSubspace for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold=1),"Threshold not identified: see ?cptSubspace for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold=PermTest))
})

test_that("Number of changepoints is correct",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,numCpts='AMOC',threshold='Manual',thresholdValue=10),"cptCovariance")
		  expect_is(cptSubspace(data2Change,subspaceDim=5,numCpts='BinSeg',threshold='Manual',thresholdValue=10),"cptCovariance")
		  expect_is(cptSubspace(data2Change,subspaceDim=5,numCpts='BinSeg',threshold='PermTest',thresholdValue=0.05),"cptCovariance")
		  expect_is(cptSubspace(data2Change,subspaceDim=5,numCpts=2,threshold='Manual',thresholdValue=10),"cptCovariance")

		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,numCpts='AMC'),"numCpts not identified: see ?cptSubspace for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,numCpts=TRUE),"numCpts not identified: see ?cptSubspace for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,numCpts=c(1,2)),"numCpts not identified: see ?cptSubspace for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,numCpts=amoc))
})

test_that("Threshold value is correct",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05),"cptCovariance")
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='Manual',thresholdValue=30),"cptCovariance")

		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,thresholdValue=-1),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='PermTest'")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=2),"When using threshold='PermTest', the thresholdValue is the significance level of the permutation test. Please enter a value between 0 and 1")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,thresholdValue="two"),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='PermTest'")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,thresholdValue=c(0.05,0.1)),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='PermTest'")
})

test_that("Minimum segment length is appropriate",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,msl=ncol(dataAMOC)+5,threshold='Manual',thresholdValue=10),"cptCovariance")

		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,msl=nrow(dataAMOC)-1),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,msl=-5),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,msl=ncol(dataAMOC)-2),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,msl='minimum'),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,msl=c(ncol(dataAMOC)+1,ncol(dataAMOC)+2)),"Minimum segment length should be a single integer between p and n/2")
})

test_that("Number of permutations is appropriate",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05,nperm=200),"cptCovariance")

		  expect_warning(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05,nperm=20),"Number of permutations is small - threshold may be unreliable")
		  expect_warning(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05,nperm=2000),"Number of permutations is large - method may take substantial time to run")
		  expect_warning(cptSubspace(dataAMOC,subspaceDim=5,threshold='Manual',nperm=100),"nperm is only used with threshold='PermTest'",fixed=TRUE)
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05,nperm=NA),"Number of permutations should be a single positive integer")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05,nperm='high'),"Number of permutations should be a single positive integer")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,threshold='PermTest',thresholdValue=0.05,nperm=c(100,200)),"Number of permutations should be a single positive integer")
})

test_that("Class argument is logical",{
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,Class=TRUE,threshold='Manual',thresholdValue=10),"cptCovariance")
		  expect_is(cptSubspace(dataAMOC,subspaceDim=5,Class=FALSE,threshold='Manual',thresholdValue=10),"integer")

		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,Class='S4'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,Class='TRUE'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptSubspace(dataAMOC,subspaceDim=5,Class=yes))
})
##}}}



