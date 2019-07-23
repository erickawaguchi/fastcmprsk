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
#' @param variance Logical: Get standard error estimates for parameter estimates via bootstrap.
#' @param var.control List of options for variance estimation.
#' @param returnDataFrame Logical: Return (ordered) data frame.
#'
#' @details Fits the 'proportional subdistribution hazards' regression model described in Fine and Gray (1999) using a novel two-way linear scan approach.
#'
#' @return Returns a list of class \code{fcrr}.
#' @importFrom survival survfit
#' @import foreach
#' @export
#' @useDynLib fastcmprsk, .registration = TRUE
#' @examples
#' library(fastcmprsk)
#'
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit <- fastCrr(ftime, fstatus, cov, variance = FALSE)
#'
#' #To parallellize variance estimation (make sure doParallel is loaded)
#' myClust <- makeCluster(2)
#' registerDoParallel(myClust)
#' fit1 <- fastCrr(ftime, fstatus, cov, variance = FALSE,
#' var.control = varianceControl(B = 100, useMultipleCores = TRUE))
#' stopCluster(myClust)
#'
#'
#'
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#'

fastCrr <- function(ftime, fstatus, X, failcode = 1, cencode = 0,
                    eps = 1E-6,
                    max.iter = 1000, getBreslowJumps = TRUE,
                    standardize = TRUE,
                    variance = TRUE,
                    var.control = varianceControl(B = 100, useMultipleCores = FALSE),
                    returnDataFrame = FALSE){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")
  if(!is.matrix(X)) X = as.matrix(X)

  # Sort time
  n <- length(ftime)
  p <- ncol(X)
  dat <- setupData(ftime, fstatus, X, cencode, failcode, standardize)

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
               as.integer(p), as.integer(n), as.double(dat$wt), as.double(denseFit[[1]] / dat$scale),
               double(sum(dat$fstatus == 1)),
               PACKAGE = "fastcmprsk")
    getBreslowJumps <- data.frame(time = unique(rev(dat$ftime[dat$fstatus == 1])),
                                  jump = as.vector(rev(unique(bjump[[8]])) * table(dat$ftime[dat$fstatus == 1], dat$fstatus[dat$fstatus == 1])))
  } #End Breslow jump

  #Calculate variance (if turned on)
  sigma <- NULL
    if(variance) {
    controls = var.control
    if (!missing(controls))
      controls[names(controls)] <- controls
    B        <- controls$B
    seed     <- controls$seed
    mcores   <- controls$mcores
    # Are we using multiple cores (parallel) or not
    if(mcores) `%mydo%` <- `%dopar%`
    else          `%mydo%` <- `%do%`

    i <- NULL  #this is only to trick R CMD check,

    set.seed(seed)
    seeds = sample.int(2^25, B, replace = FALSE)
    bsamp_beta <- numeric() #Store bootstrap values of beta here
    # %mydo% will determine whether we are using multiple cores or a single core
    bsamp_beta <- foreach(i = seeds, .combine = 'rbind', .packages = "fastcmprsk") %mydo% {
      set.seed(i)
      bsamp  <- sample(n, n, replace = TRUE) #Bootstrap sample index
      dat.bs    <- setupData(ftime[bsamp], fstatus[bsamp], X[bsamp, ], cencode, failcode, standardize)
      fit.bs <- .Call("ccd_dense", dat.bs$X, as.numeric(dat.bs$ftime), as.integer(dat.bs$fstatus), dat.bs$wt,
                      eps, as.integer(max.iter), PACKAGE = "fastcmprsk")
      tmp <- fit.bs[[1]] / dat.bs$scale
      rm(dat.bs, fit.bs)
      return(tmp)
    }
    sigma = cov(bsamp_beta) #Get variance-covariance matrix
  } #End variance option

  if(returnDataFrame) {
    df <- data.frame(ftime = ftime, fstatus = fstatus, X)
  } else {
    df <- NULL
  }

  lrt = denseFit[[2]][1] - denseFit[[2]][2] #Calculate lilkelihood ratio test
  converged <- ifelse(denseFit[[3]] < max.iter, TRUE, FALSE)
  #Results to store:
  val <- structure(list(coef = denseFit[[1]] / dat$scale,
                        var = sigma,
                        logLik = denseFit[[2]][2] / -2,
                        logLik.null = denseFit[[2]][1] / -2,
                        lrt = lrt,
                        iter = denseFit[[3]],
                        converged = converged,
                        breslowJump = getBreslowJumps,
                        uftime = unique(rev(dat$ftime[dat$fstatus == 1])),
                        isVariance = variance,
                        df = df,
                        call = sys.call()),
                   class = "fcrr")

  val
}
