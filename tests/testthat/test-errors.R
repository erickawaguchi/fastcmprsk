library("testthat")
library("fastcmprsk")

context("test-errors.R")


test_that("fastCrrp throws error for unknown penalty", {
  set.seed(10)
  ftime   <- rexp(50)
  fstatus <- sample(0:2, 50, replace = TRUE)
  cov     <- matrix(runif(250), nrow = 50)
  dimnames(cov)[[2]] <- c('x1', 'x2', 'x3', 'x4', 'x5')

  expect_that(fastCrrp(Crisk(ftime, fstatus) ~ cov, lambda = 0, penalty = "LASO"), throws_error())
})

test_that("fastCrrp throws error for negative value of lambda", {
  set.seed(10)
  ftime   <- rexp(50)
  fstatus <- sample(0:2, 50, replace = TRUE)
  cov     <- matrix(runif(250), nrow = 50)
  dimnames(cov)[[2]] <- c('x1', 'x2', 'x3', 'x4', 'x5')

  expect_that(fastCrrp(Crisk(ftime, fstatus) ~ cov, lambda = -0.1, penalty = "RIDGE"), throws_error())
})
