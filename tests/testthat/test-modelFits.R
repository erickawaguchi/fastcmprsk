library("testthat")
library("fastcmprsk")
library("cmprsk")
library("Matrix")
library("crrp")

test_that("Comapre crr with fastCrr", {
  set.seed(4291)
  ftime <- rexp(200)
  fstatus <- sample(0:2,200,replace=TRUE)
  cov <- matrix(runif(600),nrow=200)

  fit.crr    <- crr(ftime, fstatus, cov, variance = FALSE)
  fit.fast   <- fastCrr(ftime, fstatus, cov, getVariance = FALSE)
  expect_equal(as.vector(fit.crr$coef), as.vector(fit.fast$coef), tolerance = 1E-5)
})

test_that("Compare crrp with fastCrrp ", {
  set.seed(4291)
  ftime <- rexp(200)
  fstatus <- sample(0:2,200,replace=TRUE)
  cov <- matrix(runif(600),nrow=200)

  #LASSO
  fit.crrp    <- crrp(ftime, fstatus, cov, penalty = "LASSO", lambda = 0.05)
  fit.fast   <- fastCrrp(ftime, fstatus, cov, penalty = "lasso", lambda = 200 * 0.05)
  expect_equal(as.vector(fit.crrp$beta), as.vector(fit.fast$coef), tolerance = 1E-5)

  #SCAD
  fit.crrp    <- crrp(ftime, fstatus, cov, penalty = "SCAD", lambda = 0.05)
  fit.fast   <- fastCrrp(ftime, fstatus, cov, penalty = "scad", lambda = 200 * 0.05)
  expect_equal(as.vector(fit.crrp$beta), as.vector(fit.fast$coef), tolerance = 1E-5)


  #MCP
  fit.crrp    <- crrp(ftime, fstatus, cov, penalty = "MCP", lambda = 0.05)
  fit.fast   <- fastCrrp(ftime, fstatus, cov, penalty = "mcp", lambda = 200 * 0.05)
  expect_equal(as.vector(fit.crrp$beta), as.vector(fit.fast$coef), tolerance = 1E-5)


})


test_that("Compare crr with sparseCrr and fastCrr", {
  set.seed(4291)
  beta1 <- c(0.1, 0, 0.25, 0.3)
  dat <- simulateSparseFineGrayModel(nobs = 200, numCovsPerObs = 0.25, beta1 = beta1, beta2 = -beta1,
                                     parallel = FALSE)
  X <- sparseMatrix(i = dat$covariate$rowId, j = dat$covariate$covariateId,
                    x = dat$covariate$covariateValue)
  X <- as.matrix(X)

  fit.crr    <- crr(dat$outcomes$time, dat$outcomes$y, X, variance = FALSE)
  fit.fast   <- fastCrr(dat$outcomes$time, dat$outcomes$y, X, getVariance = FALSE)
  fit.sparse <- sparseCrr(dat$outcomes, dat$covariate, getVariance = FALSE)

  expect_equal(as.vector(fit.crr$coef), as.vector(fit.fast$coef),
               tolerance = 1E-5)

  expect_equal(as.vector(fit.sparse$coef), as.vector(fit.fast$coef),
               tolerance = 1E-5)
})
