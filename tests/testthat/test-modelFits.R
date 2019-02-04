library("testthat")
library("fastcmprsk")
library("cmprsk")
library("Matrix")
library("crrp")

test_that("Compare crr with fastCrr", {
  set.seed(4291)
  ftime <- rexp(200)
  fstatus <- sample(0:2,200,replace=TRUE)
  cov <- matrix(runif(600),nrow=200)

  fit.crr    <- crr(ftime, fstatus, cov, variance = FALSE)
  fit.fast   <- fastCrr(ftime, fstatus, cov, getVariance = FALSE)
  expect_equal(as.vector(fit.crr$coef), as.vector(fit.fast$coef), tolerance = 1E-5)
})

test_that("Compare crr with fastCrr w/ tied data", {
  set.seed(4291)
  ftime <- round(rexp(200) + 50, 0)
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


