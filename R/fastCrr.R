#' Fast Fine-Gray Model Estimation
#'
#' @description Estimates parameters for the proportional subdistribution hazards model using two-way linear scan approach.
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a Crisk object as returned by the \code{Crisk} function.
#' @param data a data.frame in which to interpret the variables named in the formula.
#' @param eps Numeric: algorithm stops when the relative change in any coefficient is less than \code{eps} (default is \code{1E-6})
#' @param max.iter Numeric: maximum iterations to achieve convergence (default is 1000)
#' @param getBreslowJumps Logical: Output jumps in Breslow estimator for the cumulative hazard.
#' @param standardize Logical: Standardize design matrix.
#' @param variance Logical: Get standard error estimates for parameter estimates via bootstrap.
#' @param var.control List of options for variance estimation.
#' @param returnDataFrame Logical: Return (ordered) data frame.
#'
#' @details Fits the 'proportional subdistribution hazards' regression model described in Fine and Gray (1999) using a novel two-way linear scan approach.
#' By default, the \code{Crisk} object will specify which observations are censored (0), the event of interest (1), or competing risks (2).
#'
#' @return Returns a list of class \code{fcrr}.
#' \item{coef}{the estimated regression coefficients}
#' \item{var}{estimated variance-covariance matrix via bootstrap (if \code{variance = TRUE})}
#' \item{logLik}{log-pseudo likelihood at the estimated regression coefficients}
#' \item{logLik.null}{log-pseudo likelihood when the regression coefficients are 0}
#' \item{lrt}{log-pseudo likelihood ratio test statistic for the estimated model vs. the null model.}
#' \item{iter}{iterations of coordinate descent until convergence}
#' \item{converged}{logical.}
#' \item{breslowJump}{Jumps in the Breslow baseline cumulative hazard (used by \code{predict.fcrr})}
#' \item{uftime}{vector of unique failure (event) times}
#' \item{isVariance}{logical to return if variance is chosen to be estimated}
#' \item{df}{returned ordered data frame if \code{returnDataFrame = TRUE}.}
#'
#' @importFrom survival survfit
#' @import foreach
#' @export
#' @useDynLib fastcmprsk, .registration = TRUE
#' @examples
#'
#' library(fastcmprsk)
#'
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit <- fastCrr(Crisk(ftime, fstatus) ~ cov, variance = FALSE)
#'
#' # Not run: How to set up multiple cores for boostrapping
#' # library(doParallel) #  make sure necessary packages are loaded
#' # myClust <- makeCluster(2)
#' # registerDoParallel(myClust)
#' # fit1 <- fastCrr(Crisk(ftime, fstatus) ~ cov, variance = TRUE,
#' # var.control = varianceControl(B = 100, useMultipleCores = TRUE))
#' # stopCluster(myClust)
#'
#'
#'
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#'
#'#' Kawaguchi, E.S., Shen J.I., Suchard, M. A., Li, G. (2020) Scalable Algorithms for Large Competing Risks Data, Journal of Computational and Graphical Statistics

fastCrr <- function(formula, data,
                    eps = 1E-6,
                    max.iter = 1000, getBreslowJumps = TRUE,
                    standardize = TRUE,
                    variance = TRUE,
                    var.control = varianceControl(B = 100, useMultipleCores = FALSE),
                    returnDataFrame = FALSE){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")

  # Setup formula object
  #----------
  cl <- match.call() #
  mf.all <- match.call(expand.dots = FALSE)
  m.d <- match(c("formula", "data"), names(mf.all), 0L)
  mf.d <- mf.all[c(1L, m.d)]
  mf.d$drop.unused.levels <- TRUE
  mf.d[[1L]] <- quote(stats::model.frame)
  mf.d <- eval(mf.d, parent.frame())
  outcome <- model.response(mf.d)

  # Check to see if outcome is of class Crisk
  if (!inherits(outcome, "Crisk")) stop("Outcome must be of class Crisk")
  ftime   <- as.numeric(outcome[, 1])
  fstatus <- as.numeric(outcome[, 2])

  # Design matrix
  mt.d <- attr(mf.d, "terms")

  X <- as.matrix(model.matrix(mt.d, mf.d)[, -1])
  dlabels <- labels(X)[[2]]
  #----------

  # Sort time
  n <- length(ftime)
  p <- ncol(X)
  cencode = 0; failcode = 1 #Preset
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
  bsamp_beta <- NULL
    if(variance) {
    controls = var.control
    if (!missing(controls))
      controls[names(controls)] <- controls
    B        <- controls$B
    seed     <- controls$seed
    mcores   <- controls$mcores
    extract  <- controls$extract
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
    df <- data.frame(ftime = ftime[dat$idx], fstatus = fstatus[dat$idx], X[dat$idx, ])
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
                        bootVar = bsamp_beta,
                        df = df,
                        call = sys.call()),
                   class = "fcrr")

  val
}
