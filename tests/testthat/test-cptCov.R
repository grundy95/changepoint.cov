set.seed(1)
dataAMOC <- wishartDataGeneration(n=100,p=3,tau=50)$data
dataAMOC2 <- wishartDataGeneration(n=100,p=22,tau=50)$data

test_that("Default arguments produce no errors",{
		  expect_s4_class(cptCov(X=dataAMOC,method="Ratio"),"cptCovariance")
		  expect_s4_class(cptCov(X=dataAMOC,method="CUSUM"),"cptCovariance")
})

##{{{ Error catching tests
test_that("Method argument is correct",{
		  expect_s4_class(cptCov(X=dataAMOC,method="Ratio"),"cptCovariance")
		  expect_s4_class(cptCov(X=dataAMOC,method="CUSUM"),"cptCovariance")
		  expect_s4_class(cptCov(X=dataAMOC,method="Ryan"),"cptCovariance")
		  expect_s4_class(cptCov(X=dataAMOC,method="Aue"),"cptCovariance")

		  expect_warning(cptCov(X=dataAMOC),"no method was chosen. As p<=10 the CUSUM method will be implemented")
		  expect_warning(cptCov(X=dataAMOC2),"no method was chosen. As p>10 the Ratio method will be implemented")

		  expect_error(cptCov(X=dataAMOC,method=c('subspace','Ratio')),"only one method can be implemented at once",fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method=1),"method not recognized: Please choose between 'Ratio' and 'CUSUM'",fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method='Fisher'),"method not recognized: Please choose between 'Ratio' and 'CUSUM'",fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method='Subspace'),'For subspace changepoint detection use the function cptSubspace. See ?cptSubspace for details',fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method='Grundy'),'For subspace changepoint detection use the function cptSubspace. See ?cptSubspace for details',fixed=TRUE)
		  expect_warning(cptCov(X=dataAMOC,method='Ratio',LRCov='Empirical'),'Long run covariance estimator is not used in Ratio method',fixed=TRUE)
		  expect_error(cptCov(X=dataAMOC,method=aue))
})

test_that("Data is correct format",{
		  expect_s4_class(cptCov(X=as.data.frame(dataAMOC),method="Ratio"),"cptCovariance")
})

test_that("Correct number of changepoints are returned",{
		  ans <- cptCov(X=dataAMOC,method='Ratio',numCpts=2)
		  expect_equal(length(cpts(ans))-1,2)
		  ans2 <- cptCov(X=dataAMOC,method='CUSUM',numCpts=2)
		  expect_equal(length(cpts(ans2))-1,2)
})
