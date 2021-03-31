set.seed(1)
data2Change <- wishartDataGeneration(n=200,p=5,tau=c(75,125))$data
dataNull <- matrix(rnorm(200*5),ncol=5)

data2ChangeSub <- subspaceDataGeneration(n=100,p=10,subspaceDim=3,tau=c(30,60))$data
dataNullSub <- subspaceDataGeneration(n=100,p=10,subspaceDim=3,changeSize=0)$data

test_that("Basic subspace example works",{
		  set.seed(1)
		  ans1 <- binSeg(X=data2ChangeSub, subspaceDim=3, method='Subspace',
		                 msl=10, threshold='PermTest', thresholdValue=0.05)
		  expect_type(ans1,"list")
		  expect_equal(nrow(ans1),3)

		  ans2 <- binSeg(X=data2ChangeSub,subspaceDim=3,method='Subspace',msl=10,threshold='Manual',thresholdValue=5)
		  expect_type(ans2,'list')
		  expect_equal(nrow(ans2),3)

		  ans3 <- binSeg(X=data2ChangeSub,subspaceDim=3,method='Subspace',msl=10,m=2)
		  expect_type(ans3,'list')
		  expect_equal(nrow(ans3),2)

		  ans4 <- binSeg(X=data2ChangeSub,subspaceDim=3,method='Subspace',msl=10,m=4)
		  expect_type(ans4,'list')
		  expect_equal(nrow(ans4),4)

		  set.seed(1)
		  ans5 <- binSeg(X=dataNullSub,subspaceDim=3,method='Subspace',msl=10,threshold='PermTest',thresholdValue,alpha=0.1)
		  expect_type(ans5,"list")
		  expect_equal(nrow(ans5),1)
		  expect_true(ans5$T<ans5$thresholdValue)
		  expect_false(ans5$significant)
})

test_that("Basic ratio example works",{
		  ans1 <- binSeg(X=data2Change,method='Ratio',msl=20,thresholdValue=10)
		  expect_type(ans1,"list")
		  expect_equal(nrow(ans1),3)

		  ans2 <- binSeg(X=data2Change,method='Ratio',msl=20,m=2)
		  expect_type(ans2,'list')
		  expect_equal(nrow(ans2),2)

		  ans3 <- binSeg(X=data2Change,method='Ratio',msl=20,m=4)
		  expect_type(ans3,'list')
		  expect_equal(nrow(ans3),4)

		  ans4 <- binSeg(X=dataNull,method='Ratio',thresholdValue=20,msl=20)
		  expect_type(ans4,"list")
		  expect_equal(nrow(ans4),1)
		  expect_true(ans4$T<ans4$thresholdValue)
		  expect_false(ans4$significant)
})

test_that("Basic cusum example works",{
		  ans1 <- binSeg(X=data2Change,method='CUSUM',msl=20,thresholdValue=3)
		  expect_type(ans1,"list")
		  expect_equal(nrow(ans1),3)

		  ans2 <- binSeg(X=data2Change,method='CUSUM',msl=20,m=2)
		  expect_type(ans2,'list')
		  expect_equal(nrow(ans2),2)

		  ans3 <- binSeg(X=data2Change,method='CUSUM',msl=20,m=4)
		  expect_type(ans3,'list')
		  expect_equal(nrow(ans3),4)

		  ans4 <- binSeg(X=dataNull,method='CUSUM',thresholdValue=20,msl=20)
		  expect_type(ans4,"list")
		  expect_equal(nrow(ans4),1)
		  expect_true(ans4$T<ans4$thresholdValue)
		  expect_false(ans4$significant)
})

test_that("Minimum segment length requirements trigger warnings and errors",{
		  expect_warning(binSeg(X=data2ChangeSub,subspaceDim=3,method='Subspace',msl=25,m=10),"Cannot allocate 10 changepoints due to minimum segment length restrictions",fixed=TRUE)
		  expect_error(binSeg(X=data2Change,subspaceDim=3,method='Subspace',msl=1000,m=1),"Minimum segment length should be a single integer between p and n/2")
})
