context('cptCov tests')
library(changepoint.cov)

##{{{ Data Creation
set.seed(1)
dataAMOC <- wishartDataGeneration(n=100,p=3,tau=50)$data
dataAMOC2 <- wishartDataGeneration(n=100,p=22,tau=50)$data
##}}}

##{{{ Basic Functionality
test_that("Default arguments produce no errors",{
		  expect_is(cptCov(X=dataAMOC,method="Ratio"),"cptCovariance")
		  expect_is(cptCov(X=dataAMOC,method="CUSUM"),"cptCovariance")
})
##}}}

##{{{ Error catching tests
test_that("Method argument is correct",{
		  expect_is(cptCov(X=dataAMOC,method="Ratio"),"cptCovariance")
		  expect_is(cptCov(X=dataAMOC,method="CUSUM"),"cptCovariance")
		  expect_is(cptCov(X=dataAMOC,method="Ryan"),"cptCovariance")
		  expect_is(cptCov(X=dataAMOC,method="Aue"),"cptCovariance")

		  expect_warning(cptCov(X=dataAMOC),"no method was chosen. As p<=20 the CUSUM method will be implemented")
		  expect_warning(cptCov(X=dataAMOC2),"no method was chosen. As p>20 the Ratio method will be implemented")

		  expect_error(cptCov(X=dataAMOC,method=c('subspace','Ratio')),"only one method can be implemented at once",fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method='Subspace'),'For subspace changepoint detection use the function cptSubspace. See ?cptSubspace for details',fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method='Grundy'),'For subspace changepoint detection use the function cptSubspace. See ?cptSubspace for details',fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method=aue))
})

test_that("Data is correct format",{
		  expect_is(cptCov(X=as.data.frame(dataAMOC),method="Ratio"),"cptCovariance")
})