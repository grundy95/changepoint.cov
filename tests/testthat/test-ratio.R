context('cptRatio Tests')
library(changepoint.cov)

##{{{ Data Creation
set.seed(1)
dataAMOC <- wishartDataGeneration(n=200,p=5,tau=50)$data
set.seed(1)
dataNull <- matrix(rnorm(100*3),ncol=3)
##}}}

##{{{ Basic Functionality
test_that("Default arguments produce no errors",{
		  expect_is(cptRatio(X=dataAMOC),"cptCovariance")
		  expect_is(cptRatio(X=dataNull),"cptCovariance")
})

test_that("ratioTestStat output is correct",{
		  n <- nrow(dataAMOC)
		  msl <- 2*ncol(dataAMOC)
		  ans <- ratioTestStat(dataAMOC,msl)
		  expect_is(ans,"numeric")
		  expect_equal(length(ans),n)
		  expect_equal(sum(is.na(ans)),2*msl-1)
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
		  expect_is(cptRatio(dataAMOC),"cptCovariance")
		  expect_is(cptRatio(dataAMOCdataFrame),"cptCovariance")
		  expect_error(cptRatio(dataAMOCunivariate),"Data should be a matrix")
		  expect_error(cptRatio(as.matrix(dataAMOCunivariate,ncol=1)),"Univariate time series analysis not supported")
		  expect_error(cptRatio(dataAMOCcharacter),"Data must be numeric")
		  expect_error(cptRatio(dataAMOCna),"Missing value: NA is not allowed in the data")
})

test_that("Threshold type is correct",{
		  expect_is(cptRatio(dataAMOC,threshold='Asymptotic'),"cptCovariance")
		  expect_is(cptRatio(dataAMOC,threshold='Manual'),"cptCovariance")
		  expect_error(cptRatio(dataAMOC,threshold='man'),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,threshold=1),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,threshold=Normal))
})

test_that("Number of changepoints is correct",{
		  expect_is(cptRatio(dataAMOC,numCpts='AMOC'),"cptCovariance")
		  expect_is(cptRatio(dataAMOC,numCpts='BinSeg'),"cptCovariance")
		  expect_is(cptRatio(dataAMOC,numCpts=1),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,numCpts='AMC'),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,numCpts=TRUE),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,numCpts=c(4,5)),"numCpts not identified: see?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,numCpts=amoc))
})

test_that("Threshold value is correct",{
		  expect_is(cptRatio(dataAMOC,threshold='Asymptotic'),"cptCovariance")
		  expect_is(cptRatio(dataAMOC,threshold='Manual',thresholdValue=30),"cptCovariance")

		  expect_warning(cptRatio(dataAMOC,threshold='Asymptotic',thresholdValue=0.95),"thresholdValue is not used when using threshold='Asymptotic'")
		  expect_error(cptRatio(dataAMOC,thresholdValue=-1),"thresholdValue must be a single numeric and positive if threshold='Manual'")
		  expect_error(cptRatio(dataAMOC,thresholdValue="two"),"thresholdValue must be a single numeric and positive if threshold='Manual'")
		  expect_error(cptRatio(dataAMOC,thresholdValue=c(0.05,0.1)),"thresholdValue must be a single numeric and positive if threshold='Manual'")
})

test_that("Minimum segment length is appropriate",{
		  expect_is(cptRatio(dataAMOC,msl=ncol(dataAMOC)+5),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,msl=nrow(dataAMOC)-1),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl=-5),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl=ncol(dataAMOC)-2),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl='minimum'),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl=c(ncol(dataAMOC)+1,ncol(dataAMOC)+2)),"Minimum segment length should be a single integer between p and n/2")
})

test_that("errorCheck is logical",{
		  expect_is(cptRatio(dataAMOC,errorCheck=TRUE),"cptCovariance")
		  expect_is(cptRatio(dataAMOC,errorCheck=FALSE),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,errorCheck='FALSE'),"errorCheck should be logical")
})

test_that("Class argument is logical",{
		  expect_is(cptRatio(dataAMOC,Class=TRUE),"cptCovariance")
		  expect_is(cptRatio(dataAMOC,Class=FALSE),"integer")

		  expect_error(cptRatio(dataAMOC,Class='S4'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptRatio(dataAMOC,Class='TRUE'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptRatio(dataAMOC,Class=yes))
})
##}}}

