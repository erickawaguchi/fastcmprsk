#' Fine-Gray Model Estimation
#'
#' @description Estimates parameters for the proportional subdistribution hazards model.
#' Penalties currently include LASSO and ridge regression. User-specificed weights can be assigned
#' to the penalty for each coefficient (e.g. implementing adaptive LASSO).
#'
#' @param ftime A vector of event/censoring times.
#' @param fstatus A vector with unique code for each event type and a separate code for censored observations.
#' @param X A matrix of fixed covariates (nobs x ncovs)
#' @param failcode Integer: code of \code{fstatus} that event type of interest (default is 1)
#' @param cencode Integer: code of \code{fstatus} that denotes censored observations (default is 0)
#' @param eps Numeric: algorithm stops when the relative change in any coefficient is less than \code{eps} (default is \code{1E-6})
#' @param max.iter Numeric: maximum iterations to achieve convergence (default is 1000)
#' @param getBreslowJumps Logical: Output jumps in Breslow estimator for the cumulative hazard
#' @param getVariance Logical: Get standard error estimates for parameter estimates via bootstrap or quadratic approximation.
#' @param returnDataFrame Logical: Return data frame.
#'
#' @details The \code{pshBAR} function penalizes the log-partial likelihood of the proportional subdistribution hazards model
#' from Fine and Gray (1999) with the Broken Adaptive Ridge (BAR) penalty. A cyclic coordinate descent algorithm is used for implementation.
#' For stability, the covariate matrix \code{X} is standardized prior to implementation.
#'
#' Special cases: Fixing \code{xi} and \code{lambda} to 0 results in the standard competing risk regression using \code{crr}.
#' Fixing \code{lambda} to 0 and specifying \code{xi} will result in a ridge regression solution.
#' @return Returns a list of class \code{pshBAR}.
#'
#' @importFrom survival survfit
#' @import doParallel
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

fastCrr <- function(ftime, fstatus, X, failcode = 1, cencode = 0,
                    eps = 1E-6,
                    max.iter = 1000, getBreslowJumps = TRUE,
                    standardize = TRUE,
                    startingCoefficients = NULL,
                    getVariance = TRUE,
                    varianceControl = varianceControl(...),
                    returnDataFrame = FALSE){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")

  # Sort time
  n <- length(ftime)
  p <- ncol(X)
  dat <- setupData(ftime, fstatus, X, cencode, failcode, standardize)

  scale = dat$scale
  ## Fit the PSH Ridge Model here w/ tuning parameter xi
  denseFit   <- .Call("ccd_dense", dat$X, as.numeric(dat$ftime), as.integer(dat$fstatus), dat$wt,
                      eps, as.integer(max.iter), PACKAGE = "fastcmprsk")

  if (denseFit[[3]] == max.iter) {
    warning("Maximum number of iterations reached. Estimates may not be stable")
  }

  #Calculate Breslow Baseline
  if(getBreslowJumps) {
    bjump = .C("getBreslowJumps", as.double(dat$ftime), as.integer(dat$fstatus),
               as.double(sweep(sweep(dat$X, 2, dat$scale, "*"), 2, dat$center, `+`)),
               as.integer(p), as.integer(n), as.double(dat$wt), as.double(denseFit[[1]] / scale),
               double(sum(dat$fstatus == 1)),
               PACKAGE = "fastcmprsk")
    getBreslowJumps <- data.frame(time = unique(rev(dat$ftime[dat$fstatus == 1])),
                                  jump = as.vector(rev(unique(bjump[[8]])) * table(dat$ftime[dat$fstatus == 1], dat$fstatus[dat$fstatus == 1])))
  } #End Breslow jump

  #Calculate variance (if turned on)
  method = var = NA #Set se & method to NA, will update if getVariance
  if(getVariance) {
    controls <- varianceControl
    if (!missing(controls))
      controls[names(controls)] <- controls
    B        <- controls$B
    parallel <- controls$parallel
    ncores   <- controls$ncores
    seed     <- controls$seed
    method   <- controls$method

    if(method == "approx") {
      var = sweep(solve(t(apply(dat$X, 2, function(x) x *  denseFit[[6]])) %*% dat$X), 2, dat$scale^2, FUN = "/")
    } else if (method == "bootstrap") {
      #Run if parallel
      if(parallel) {
        set.seed(seed)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        bsamp_beta <- foreach(i = 1:B, .combine = rbind) %do% {
          bsamp  <- sample(n, n, replace = TRUE) #Bootstrap sample index
          dat <- setupData(ftime[bsamp], fstatus[bsamp], X[bsamp, ], cencode, failcode, standardize)
          fit.bs   <- .Call("ccd_dense", dat$X, as.numeric(dat$ftime), as.integer(dat$fstatus), dat$wt,
                            eps, as.integer(max.iter), PACKAGE = "fastcmprsk")
          fit.bs[[1]] / dat$scale
        }
        stopCluster(cl)
      } else {
        set.seed(seed)
        bsamp_beta <- matrix(NA, nrow = B, ncol = p)
        for(i in 1:B) {
          bsamp  <- sample(n, n, replace = TRUE) #Bootstrap sample index
          dat.bs    <- setupData(ftime[bsamp], fstatus[bsamp], X[bsamp, ], cencode, failcode, standardize)
          fit.bs <- .Call("ccd_dense", dat.bs$X, as.numeric(dat.bs$ftime), as.integer(dat.bs$fstatus), dat.bs$wt,
                              eps, as.integer(max.iter), PACKAGE = "fastcmprsk")
          bsamp_beta[i, ] <- fit.bs[[1]] / dat.bs$scale
        }
      }
      #Calculate standard error
      var = cov(bsamp_beta)
    } #End bootstrap option
  } #End variance


  if(returnDataFrame) {
    df <- data.frame(ftime = ftime, fstatus = fstatus, X)
  } else {
    df <- NULL
  }

  #Results to store:
  val <- structure(list(coef = denseFit[[1]] / scale,
                        var = var,
                        logLik = denseFit[[2]] / -2,
                        iter = denseFit[[3]],
                        breslowJump = getBreslowJumps,
                        uftime = unique(rev(dat$ftime[dat$fstatus == 1])),
                        method = method,
                        df = df,
                        call = sys.call()),
                   class = "fcrr")

  val
}
