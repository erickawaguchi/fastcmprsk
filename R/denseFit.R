#' Fast Fine-Gray Model Estimation
#'
#' @description Estimates parameters for the proportional subdistribution hazards model using two-way linear scan approach.
#'
#' @param ftime A vector of event/censoring times.
#' @param fstatus A vector with unique code for each event type and a separate code for censored observations.
#' @param X A matrix of fixed covariates (nobs x ncovs)
#' @param failcode Integer: code of \code{fstatus} that event type of interest (default is 1)
#' @param cencode Integer: code of \code{fstatus} that denotes censored observations (default is 0)
#' @param eps Numeric: algorithm stops when the relative change in any coefficient is less than \code{eps} (default is \code{1E-6})
#' @param max.iter Numeric: maximum iterations to achieve convergence (default is 1000)
#' @param getBreslowJumps Logical: Output jumps in Breslow estimator for the cumulative hazard.
#' @param standardize Logical: Standardize design matrix.
#' @param getVariance Logical: Get standard error estimates for parameter estimates via bootstrap.
#' @param var.control List of options for variance estimation.
#' @param returnDataFrame Logical: Return (ordered) data frame.
#'
#' @details Fits the 'proportional subdistribution hazards' regression model described in Fine and Gray (1999) using a novel two-way linear scan approach.
#'
#' @return Returns a list of class \code{fcrr} with the following components:
#' @param $coef Fine-Gray regression coefficient estimates
#' @param $var Variance-covariance estimates via bootstrap (if \code{getVariance = TRUE})
#' @param $loglik log pseudo-likelihood evaluated at \code{$coef}.
#' @param $iter Number of iterations it took for convergence
#' @param $breslowJump Jumps in the Breslow-type estimate of the underlying sub-distribution cumulative hazard.
#' @param $uftime Vector of unique \code{failcode} event times.
#' @param $df Returned data frame that is ordered by decreasing event time. (If \code{returnDataFrame = TRUE}).
#' @importFrom survival survfit
#' @import doParallel
#' @export
#' @useDynLib fastcmprsk
#' @examples
#' library(cmprsk)
#' library(fastcmprsk)
#'
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit1 <- fastCrr(ftime, fstatus, cov, getVariance = FALSE)
#' fit2 <- crr(ftime, fstatus, cov)
#' max(abs(fit1$coef - fit2$coef))
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#'

fastCrr <- function(ftime, fstatus, X, failcode = 1, cencode = 0,
                    eps = 1E-6,
                    max.iter = 1000, getBreslowJumps = TRUE,
                    standardize = TRUE,
                    getVariance = TRUE,
                    var.control = varianceControl(B = 100),
                    returnDataFrame = FALSE){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")

  # Sort time
  n <- length(ftime)
  p <- ncol(X)
  dat <- setupData(ftime, fstatus, X, cencode, failcode, standardize)

  scale = dat$scale

  #Fit model here
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
  var = rep(NA, p) #Set se & method to NA, will update if getVariance
  if(getVariance) {
    controls = var.control
    if (!missing(controls))
      controls[names(controls)] <- controls
    B        <- controls$B
    parallel <- controls$parallel
    ncores   <- controls$ncores
    seed     <- controls$seed

      if(parallel) {
        set.seed(seed)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        bsamp_beta <- foreach(i = 1:B, .combine = rbind) %dopar% {
          bsamp  <- sample(n, n, replace = TRUE) #Bootstrap sample index
          dat <- setupData(ftime[bsamp], fstatus[bsamp], X[bsamp, ], cencode, failcode, standardize)
          fit.bs   <- .Call("ccd_dense", dat$X, as.numeric(dat$ftime), as.integer(dat$fstatus), dat$wt,
                            eps, as.integer(max.iter), PACKAGE = "fastcmprsk")
          if (fit.bs[[3]] == max.iter) {
            warning(paste0("Maximum number of iterations reached for ", i, "th bootstrap sample. Estimates may not be stable"))
          }
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
          if (fit.bs[[3]] == max.iter) {
            warning(paste0("Maximum number of iterations reached for ", i, "th bootstrap sample. Estimates may not be stable"))
          }
          bsamp_beta[i, ] <- fit.bs[[1]] / dat.bs$scale
        }
      }
      #Calculate standard error
      var = cov(bsamp_beta)
    } #End variance option

  if(returnDataFrame) {
    df <- data.frame(ftime = ftime, fstatus = fstatus, X)
  } else {
    df <- NULL
  }

  converged <- ifelse(denseFit[[3]] < max.iter, TRUE, FALSE)
  #Results to store:
  val <- structure(list(coef = denseFit[[1]] / scale,
                        var = var,
                        logLik = denseFit[[2]] / -2,
                        iter = denseFit[[3]],
                        converged = converged,
                        breslowJump = getBreslowJumps,
                        uftime = unique(rev(dat$ftime[dat$fstatus == 1])),
                        getVariance = getVariance,
                        df = df,
                        call = sys.call()),
                   class = "fcrr")

  val
}
