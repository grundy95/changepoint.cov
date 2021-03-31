set.seed(1)
dataAMOCsub <- subspaceDataGeneration(n=100,p=20,subspaceDim=5,tau=50,changeSize=0.5*sqrt(5))$data
ansSubspace <- cptSubspace(X=dataAMOCsub,subspaceDim=5,msl=20)
ansSubspaceMan <- cptSubspace(X=dataAMOCsub,numCpts=1,subspaceDim=5,msl=20)

set.seed(1)
dataAMOC <- wishartDataGeneration(n=100,p=3,tau=50)$data
ansCUSUM <- cptCUSUM(X=dataAMOC)
ansRatio <- cptRatio(X=dataAMOC)
ansCovRatio <- cptCov(X=dataAMOC,method="Ratio")
ansCovCUSUM <- cptCov(X=dataAMOC,method="CUSUM")

test_that("Object is of correct class",{
	  expect_s4_class(ansSubspace,"cptCovariance")
	  expect_s4_class(ansCUSUM,"cptCovariance")
	  expect_s4_class(ansRatio,"cptCovariance")
	  expect_s4_class(ansCovCUSUM,"cptCovariance")
	  expect_s4_class(ansCovRatio,"cptCovariance")
})

test_that("Methods show and summary are working",{
  expect_snapshot(summary(ansSubspace))
  expect_snapshot(summary(ansCUSUM))
  expect_snapshot(summary(ansRatio))
  expect_snapshot(summary(ansCovCUSUM))
  expect_snapshot(summary(ansCovRatio))

  expect_type(invisible(capture.output(show(ansSubspace))),'character')
  expect_type(invisible(capture.output(show(ansCUSUM))),'character')
  expect_type(invisible(capture.output(show(ansRatio))),'character')
  expect_type(invisible(capture.output(show(ansCovCUSUM))),'character')
  expect_type(invisible(capture.output(show(ansCovRatio))),'character')
})

test_that("Plotting is working",{
		  expect_s3_class(plot(ansSubspace),'ggplot')
		  expect_s3_class(plot(ansRatio),'ggplot')
		  expect_s3_class(plot(ansCUSUM),'ggplot')
		  expect_s3_class(plot(ansCovCUSUM),'ggplot')
		  expect_s3_class(plot(ansCovRatio),'ggplot')
})

test_that("Slots are correct",{
		  expect_type(ansSubspace@dataset,"double")
		  expect_type(ansSubspace@cpts,"integer")
		  expect_type(ansSubspace@method,"character")
		  expect_type(ansSubspace@msl,"integer")
		  expect_type(ansSubspace@numCpts,"character")
		  expect_type(ansSubspace@threshold,"character")
		  expect_type(ansSubspace@thresholdValue,"double")
		  expect_type(ansSubspace@cptsSig,"list")
		  expect_type(ansSubspace@subspaceDim,"double")
		  expect_type(ansSubspace@nperm,"double")
		  expect_type(ansSubspace@LRCov,"character")
		  expect_type(ansSubspace@date,"character")
		  expect_type(ansSubspace@version,"character")

		  expect_type(ansCUSUM@dataset,"double")
		  expect_type(ansCUSUM@cpts,"integer")
		  expect_type(ansCUSUM@method,"character")
		  expect_type(ansCUSUM@msl,"integer")
		  expect_type(ansCUSUM@numCpts,"character")
		  expect_type(ansCUSUM@threshold,"character")
		  expect_type(ansCUSUM@thresholdValue,"double")
		  expect_type(ansCUSUM@cptsSig,"list")
		  expect_type(ansCUSUM@subspaceDim,"double")
		  expect_type(ansCUSUM@nperm,"double")
		  expect_type(ansCUSUM@LRCov,"character")
		  expect_type(ansCUSUM@date,"character")
		  expect_type(ansCUSUM@version,"character")

		  expect_type(ansRatio@dataset,"double")
		  expect_type(ansRatio@cpts,"integer")
		  expect_type(ansRatio@method,"character")
		  expect_type(ansRatio@msl,"integer")
		  expect_type(ansRatio@numCpts,"character")
		  expect_type(ansRatio@threshold,"character")
		  expect_type(ansRatio@thresholdValue,"double")
		  expect_type(ansRatio@cptsSig,"list")
		  expect_type(ansRatio@subspaceDim,"double")
		  expect_type(ansRatio@nperm,"double")
		  expect_type(ansRatio@LRCov,"character")
		  expect_type(ansRatio@date,"character")
		  expect_type(ansRatio@version,"character")

		  expect_type(ansCovCUSUM@dataset,"double")
		  expect_type(ansCovCUSUM@cpts,"integer")
		  expect_type(ansCovCUSUM@method,"character")
		  expect_type(ansCovCUSUM@msl,"integer")
		  expect_type(ansCovCUSUM@numCpts,"character")
		  expect_type(ansCovCUSUM@threshold,"character")
		  expect_type(ansCovCUSUM@thresholdValue,"double")
		  expect_type(ansCovCUSUM@cptsSig,"list")
		  expect_type(ansCovCUSUM@subspaceDim,"double")
		  expect_type(ansCovCUSUM@nperm,"double")
		  expect_type(ansCovCUSUM@LRCov,"character")
		  expect_type(ansCovCUSUM@date,"character")
		  expect_type(ansCovCUSUM@version,"character")

		  expect_type(ansCovRatio@dataset,"double")
		  expect_type(ansCovRatio@cpts,"integer")
		  expect_type(ansCovRatio@method,"character")
		  expect_type(ansCovRatio@msl,"integer")
		  expect_type(ansCovRatio@numCpts,"character")
		  expect_type(ansCovRatio@threshold,"character")
		  expect_type(ansCovRatio@thresholdValue,"double")
		  expect_type(ansCovRatio@cptsSig,"list")
		  expect_type(ansCovRatio@subspaceDim,"double")
		  expect_type(ansCovRatio@nperm,"double")
		  expect_type(ansCovRatio@LRCov,"character")
		  expect_type(ansCovRatio@date,"character")
		  expect_type(ansCovRatio@version,"character")
})

test_that("Slot retrival and replacement functions are working",{
  ansSubspace@dataset <- matrix(0,100*20,ncol=20)
		  expect_type(dataset(ansSubspace),"double")
		  ansSubspace@cpts <- 50L
		  expect_type(cpts(ansSubspace),"integer")
		  ansSubspace@method <- 'Subspace'
		  expect_type(method(ansSubspace),"character")
		  ansSubspace@msl <- 10L
		  expect_type(msl(ansSubspace),"integer")
		  ansSubspace@numCpts <- 'BinSeg'
		  expect_type(numCpts(ansSubspace),"character")
		  ansSubspace@threshold <- 'PermTest'
		  expect_type(threshold(ansSubspace),"character")
		  ansSubspace@thresholdValue <- 10
		  expect_type(thresholdValue(ansSubspace),"double")
		  ansSubspace@cptsSig <- data.frame('Cpts'=100)
		  expect_type(cptsSig(ansSubspace),"list")
		  ansSubspace@subspaceDim <- 5
		  expect_type(subspaceDim(ansSubspace),"double")
		  ansSubspace@nperm <- 100
		  expect_type(nperm(ansSubspace),"double")
		  expect_type(covEst(ansSubspace),"list")
		  expect_type(subspaceEst(ansSubspace),"list")

		  ansCUSUM@LRCov <- 'Bartlett'
		  expect_type(LRCov(ansCUSUM),'character')

		  ansSubspaceMan@numCpts = 'Manual'
		  expect_type(numCpts(ansSubspaceMan),'character')

		  expect_error(subspaceDim(ansRatio),"subspaceDim is only a valid slot for method='Subspace'",fixed=TRUE)
		  expect_error(subspaceEst(ansRatio),"Subspace estimation only possible for method='Subspace'",fixed=TRUE)
		  expect_error(nperm(ansRatio),"nperm is only a valid slot when using the permutation test within method='Subspace'",fixed=TRUE)
		  expect_error(nperm(ansSubspaceMan),"nperm is only a valid slot when using the permutation test within method='Subspace'",fixed=TRUE)
		  expect_error(LRCov(ansSubspace),"LRCov is only a valid slot for method='CUSUM'")
})



