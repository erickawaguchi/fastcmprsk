#' Control for Variance Calculation
#'
#' @description Fits broken adaptive ridge regression for competing risks regression.
#' This package allows for ridge and broken adaptive ridge penalties.
#'
#' @param ftime A vector of event/censoring times.
#' @param fstatus A vector with unique code for each event type and a separate code for censored observations.
#' @param X A matrix of fixed covariates (nobs x ncovs)
#' @param failcode Integer: code of \code{fstatus} that event type of interest (default is 1)
#' @param cencode Integer: code of \code{fstatus} that denotes censored observations (default is 0)
#' @param lambda Numeric: BAR tuning parameter value
#' @param xi Numeric: tuning parameter for initial ridge regression
#' @param delta Numeric: change from 2 in ridge norm dimension
#' @param eps Numeric: algorithm stops when the relative change in any coefficient is less than \code{eps} (default is \code{1E-6})
#' @param tol Numeric: absolute threshold at which to force coefficients to 0 (default is \code{1E-6})
#' @param lam.min Numeric: smallest value of lambda if performing grid search
#' @param nlambda Numeric: number of \code{lambda} values if performing grid search  (default is 25)
#' @param log Logical: Whether or not the grid search is log10 spaced (default is \code{TRUE})
#' @param max.iter Numeric: maximum iterations to achieve convergence (default is 1000)
#'
#' @details The \code{pshBAR} function penalizes the log-partial likelihood of the proportional subdistribution hazards model
#' from Fine and Gray (1999) with the Broken Adaptive Ridge (BAR) penalty. A cyclic coordinate descent algorithm is used for implementation.
#' For stability, the covariate matrix \code{X} is standardized prior to implementation.
#'
#' Special cases: Fixing \code{xi} and \code{lambda} to 0 results in the standard competing risk regression using \code{crr}.
#' Fixing \code{lambda} to 0 and specifying \code{xi} will result in a ridge regression solution.
#' @return Returns a list of class \code{pshBAR}.
#'
#' @import survival
#' @export
#' @useDynLib fastcmprsk
#' @examples
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit <- pshBAR(ftime, fstatus, cov, lambda = log(5) / 2, xi = log(5))
#' fit$coef
#' @references
#' Breheny, P. and Huang, J. (2011) Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection. \emph{Ann. Appl. Statist.}, 5: 232-253.
#'
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.

varianceControl <- function(B = 100L, parallel = FALSE, ncores = 1, seed = 1991L,
                            method = "bootstrap")
{

  if(parallel) {
    if (ncores < 0L || ncores > detectCores() || is.null(ncores)) {
      warning("The value of 'ncores' was either 0, larger than the number of cores available of NULL. Set to 1")
      ncores <- 1L
    }
  }

  if (B <= 0) {
    warning("The value of 'B' must be a non-negative integer. Set to 100")
    B <- 100L
  }

  if(seed <= 0) {
    warnings("The value of 'seed' must be a non-negative integer. Set to 1991")
    seed <- 1991L
  }

  if(!(method %in% c("bootstrap", "approx"))) {
    stop("The value of 'method' must be either 'bootstrap' or 'approx'.")
  }
  list(B = B, parallel = parallel, ncores = ncores, seed = seed, method = method)
}
