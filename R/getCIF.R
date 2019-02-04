#' Cumulative Incidence Function Estimation
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
#' @import survival dynpred
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
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.

getCIF <- function(fit, cov, getBootstrapVariance = TRUE, B = 100,
                   type = "none",
                   alpha = 0.05, seednum = 1991, eps = 1E-6,
                   tL = NULL, tU = NULL,
                   max.iter = 1000){

  ## Error checking
  if(class(fit) != "fcrr") {
    stop("Object 'fit' must be of class fcrr")
  }

  if(is.null(fit$breslowJump)) {
    stop("Breslow jumps were not calculated. Please re-run model with 'getBreslowJumps = TRUE'")
  }

  if(is.null(fit$df)) {
    stop("Bresow jumps were not calculated. Please re-run model with 'getBreslowJumps = TRUE'")
  }

  if(!(type %in% c("none", "bands", "point"))) {
    type = "none"
    warning("type is incorrectly specified. Valid options are 'bands', 'point', 'none'.
            Set to 'none'")
  }

  if(alpha <= 0 | alpha >= 1) {
    alpha = 0.05
    warning("alpha is incorrectly specified. Set to 0.05")
  }

  if(B < 0) {
    B = 100
    warning("B is incorrectly specified. Set to 100")
  }

  if(is.null(tL)) tL <- min(fit$uftime)
  if(is.null(tU)) tU <- max(fit$uftime)

  if(tL <= 0 | tL >= max(fit$uftime)) {
    tL <- min(fit$uftime)
    warning("tL is incorrectly specified (can not be nonpositive or larger than largest observed event time.
            Set to smallest observed event time")
  }

  if(tU <= 0 | tU <= min(fit$uftime)) {
    tL <- min(fit$uftime)
    warning("tU is incorrectly specified (can not be nonpositive or smaller than smallest observed event time.
            Set to largest observed event time")
  }

  min.idx = min(which(fit$uftime >= tL))
  max.idx = max(which(fit$uftime <= tU))

  if (length(fit$coef) == length(cov)) {
    CIF.hat <- cumsum(exp(sum(cov * fit$coef)) * fit$breslowJump[, 2]) #This is cumulative hazard
    CIF.hat <- 1 - exp(-CIF.hat)
  } else {
    stop("Parameter dimension of 'cov' does not match dimension of '$coef' from object.")
  }

  #Get SD via bootstrap
  set.seed(seednum)
  if(getBootstrapVariance) {
    CIF.boot <- matrix(NA, nrow = B, ncol = length(fit$uftime))
    colnames(CIF.boot) <- round(fit$uftime, 3)
    ftime <- fit$df$ftime
    fstatus <- fit$df$fstatus
    n <- length(ftime)
    X <- as.matrix(fit$df[, -(1:2)])
    for(i in 1:B) {
      bsamp  <- sample(n, n, replace = TRUE) #Bootstrap sample index
      fit.bs <- fastCrr(ftime[bsamp], fstatus[bsamp], X[bsamp, ], getVariance = FALSE)
      CIF.bs <- 1 - exp(-cumsum(exp(sum(cov * fit.bs$coef)) * fit.bs$breslowJump[, 2]))
      CIF.boot[i, ] <- evalstep(fit.bs$breslowJump$time,
                               stepf = CIF.bs,
                               subst = 1E-16,
                               newtime = fit$uftime)
    }
    rm(CIF.bs)
  } #End bootstrap variance

  #Variance Stabalization: f(x) = log(-log(x))
  CIF.hat  <- log(-log(CIF.hat))
  CIF.boot <- log(-log(CIF.boot))
  CIF.sd <- apply(CIF.boot, 2, sd)
  if(type == "bands") {
    #If interval type is confidence band.
    #Find Pr(sup_[tL, tU] |Fhat - F| / sd(F) <= C) = 1 - alpha / 2
    sup    <- apply(CIF.boot, 1, function(x) max((abs(x - CIF.hat) / CIF.sd)[min.idx:max.idx]))
    z.stat <- quantile(sup, 1 - alpha / 2) #Find bootstrap quantile of sup|Fhat - F|
    llim   <- CIF.hat + z.stat * CIF.sd
    ulim   <- CIF.hat - z.stat * CIF.sd
    res  <- data.frame(ftime = fit$uftime, CIF = exp(-exp(CIF.hat)), lower = exp(-exp(llim)), upper = exp(-exp(ulim)))
  } else if (type == "point") {
    #If interval type if pointwise
    llim   <- CIF.hat + qnorm(1 - alpha / 2) * CIF.sd
    ulim   <- CIF.hat - qnorm(1 - alpha / 2) * CIF.sd
    res  <- data.frame(ftime = fit$uftime, CIF = exp(-exp(CIF.hat)), lower = exp(-exp(llim)), upper = exp(-exp(ulim)))
  } else {
    res  <- data.frame(ftime = fit$uftime, CIF = CIF.hat)
  }
  #Subset corresponding to tL and tU
  res <- subset(res, res$ftime >= tL & res$ftime <= tU)
  return(res)
}
