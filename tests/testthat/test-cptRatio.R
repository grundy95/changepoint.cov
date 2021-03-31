set.seed(1)
dataAMOC <- wishartDataGeneration(n=200,p=5,tau=50)$data
dataNull <- matrix(rnorm(100*3),ncol=3)
data2Change <- wishartDataGeneration(n=300,p=30,tau=c(100,200))$data
test_that("Default arguments produce no errors",{
		  expect_s4_class(cptRatio(X=dataAMOC),"cptCovariance")
		  expect_s4_class(cptRatio(X=dataNull),"cptCovariance")
})

test_that("ratioTestStat output is correct",{
		  n <- nrow(dataAMOC)
		  msl <- 2*ncol(dataAMOC)
		  ans <- ratioTestStat(dataAMOC,msl)
		  expect_type(ans,"double")
		  expect_equal(length(ans),n)
		  expect_equal(sum(is.na(ans)),2*msl-1)
})

test_that("ratioTestStat works for small msl",{
		  n <- nrow(dataAMOC)
		  msl <- ncol(dataAMOC)
		  expect_warning(ratioTestStat(dataAMOC,msl),'msl is too short making the test statistic incalculable for certain boundary changepoint locations. Method has continued with the exclusion of the incalculable potential changepoint locations.')
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
		  dataConstant <- matrix(0,ncol=5,nrow=100)
		  expect_s4_class(cptRatio(dataAMOC),"cptCovariance")
		  expect_s4_class(cptRatio(dataAMOCdataFrame),"cptCovariance")
		  expect_error(cptRatio(dataAMOCunivariate),"Data should be a matrix")
		  expect_error(cptRatio(as.matrix(dataAMOCunivariate,ncol=1)),"Univariate time series analysis not supported")
		  expect_error(cptRatio(dataAMOCcharacter),"Data must be numeric")
		  expect_error(cptRatio(dataAMOCna),"Missing value: NA is not allowed in the data")
		  expect_error(cptCUSUM(dataConstant),"Sample covariance of whole data is singular. This is probably due to constant data",fixed=TRUE)
})

test_that("Threshold type is correct",{
		  expect_s4_class(cptRatio(dataAMOC,threshold='Asymptotic'),"cptCovariance")
		  expect_s4_class(cptRatio(dataAMOC,threshold='Manual'),"cptCovariance")
		  expect_error(cptRatio(dataAMOC,threshold='man'),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,threshold=1),"Threshold not identified: see ?cptCov for valid entries to threshold. NOTE thresholds should be character strings",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,threshold=Normal))
})

test_that("Number of changepoints is correct",{
		  expect_s4_class(cptRatio(dataAMOC,numCpts='AMOC'),"cptCovariance")
		  expect_s4_class(cptRatio(data2Change,numCpts='BinSeg'),"cptCovariance")
		  expect_s4_class(cptRatio(data2Change,numCpts=2),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,numCpts='AMC'),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,numCpts=TRUE),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,numCpts=c(4,5)),"numCpts not identified: see ?cptCov for valid entries to numCpts",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,numCpts=amoc))
})

test_that("Threshold value is correct",{
		  expect_s4_class(cptRatio(dataAMOC,threshold='Asymptotic'),"cptCovariance")
		  expect_s4_class(cptRatio(dataAMOC,threshold='Manual',thresholdValue=30),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,thresholdValue=-1),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,thresholdValue="two"),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'",fixed=TRUE)
		  expect_error(cptRatio(dataAMOC,thresholdValue=c(0.05,0.1)),"thresholdValue must be a single numeric and positive if threshold='Manual' or between 0 and 1 if threshold='Asymptotic'",fixed=TRUE)
})

test_that("Minimum segment length is appropriate",{
		  expect_s4_class(cptRatio(dataAMOC,msl=ncol(dataAMOC)+5),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,msl=nrow(dataAMOC)-1),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl=-5),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl=ncol(dataAMOC)-2),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl='minimum'),"Minimum segment length should be a single integer between p and n/2")
		  expect_error(cptRatio(dataAMOC,msl=c(ncol(dataAMOC)+1,ncol(dataAMOC)+2)),"Minimum segment length should be a single integer between p and n/2")
})

test_that("errorCheck is logical",{
		  expect_s4_class(cptRatio(dataAMOC,errorCheck=TRUE),"cptCovariance")
		  expect_s4_class(cptRatio(dataAMOC,errorCheck=FALSE),"cptCovariance")

		  expect_error(cptRatio(dataAMOC,errorCheck='FALSE'),"errorCheck should be logical")
})

test_that("Class argument is logical",{
		  expect_s4_class(cptRatio(dataAMOC,Class=TRUE),"cptCovariance")
		  expect_type(cptRatio(dataAMOC,Class=FALSE),"integer")

		  expect_error(cptRatio(dataAMOC,Class='S4'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptRatio(dataAMOC,Class='TRUE'),"Class should be logical, TRUE or FALSE")
		  expect_error(cptRatio(dataAMOC,Class=yes))
})
##}}}

