library("testthat")
library("fastcmprsk")

context("test-ENET.R")


test_that("Test Elastic Net (LASSO)", {
  set.seed(4291)
  ftime <- rexp(200)
  fstatus <- sample(0:2,200,replace=TRUE)
  cov <- matrix(runif(600),nrow=200)

  fit.lasso <- fastCrrp(Crisk(ftime, fstatus) ~ cov, penalty = "LASSO", lambda = 0.01)
  fit.enet  <- fastCrrp(Crisk(ftime, fstatus) ~ cov, penalty = "ENET", alpha = 1,  lambda = 0.01)
  expect_equal(as.vector(fit.lasso$coef), as.vector(fit.enet$coef), tolerance = 1E-8)

  fit.ridge <- fastCrrp(Crisk(ftime, fstatus) ~ cov, penalty = "RIDGE", lambda = 0.05)
  fit.enet  <- fastCrrp(Crisk(ftime, fstatus) ~ cov, penalty = "ENET", alpha = 0,  lambda = 0.05)
  expect_equal(as.vector(fit.ridge$coef), as.vector(fit.enet$coef), tolerance = 1E-8)
})
