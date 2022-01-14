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

test_that("Crisk correctly categorizes event types", {
  ftime <- 1:10
  fstatus <- c(0, 1, 0, 2, 2, 1, 1, 1, 0, 0)

  test1 <- Crisk(ftime, fstatus, cencode = 2, failcode = 0)
  expect_equivalent(test1[, 2], c(1, 2, 1, 0, 0, 2, 2, 2, 1, 1))

  test1 <- Crisk(ftime, fstatus, cencode = 1, failcode = 2)
  expect_equivalent(test1[, 2], c(2, 0, 2, 1, 1, 0, 0, 0, 2, 2))
})

