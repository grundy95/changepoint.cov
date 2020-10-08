context('data Generation Tests')
library(changepoint.cov)

test_that("subspace data generation is working",{
		  expect_is(subspaceDataGeneration(n=100,p=10,q=7,tau=c(30,60),changeSize=0.5*sqrt(3)),"list")
		  expect_is(subspaceDataGeneration(n=100,p=10,q=1,tau=c(30,60),changeSize=0.5*sqrt(3)),"list")

		  expect_error(subspaceDataGeneration(n=100,p=10,q=3,tau=c(25,50,75),changeSize=c(0.5*sqrt(3),0.4*sqrt(3))),"Number of changepoints does not match number of change sizes",fixed=TRUE)
		  expect_error(subspaceDataGeneration(n=100,p=10,q=3,tau=50,changeSize=1.1*sqrt(3)),"Change sizes must be between 0 and sqrt(q)",fixed=TRUE)
})

test_that("Wishart data generation is working",{
		  Sigma <- list("sigma1"=diag(rep(1,3)),"sigma2"=diag(rep(2,3)))
		  expect_error(wishartDataGeneration(n=100,p=3,tau=c(25,50,75),Sigma=Sigma),"Need 1 covariance matrix per segment")
		  expect_is(wishartDataGeneration(n=100,p=3),"list")
})




