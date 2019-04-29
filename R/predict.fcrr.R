#' Cumulative Incidence Function Estimation
#'
#' @description Predicts cumulative incidence function from a \code{fcrr} object.
#'
#' @param fit Output from \code{fcrr} object.
#' @param cov A set of covariate values.
#' @param getBootstrapVariance Logical: Calculate variance for CIF via bootstrap.
#' @param B Number of bootstrap samples for variance estimation.
#' @param type Confidence intervals or confidence bands.
#' @param alpha Significance level to compute intervals or bands
#' @param seed Seed number of bootstrap sampling.
#' @param tL Lower time for band estimation.
#' @param tU Upper time for band estimation.
#' @param ... additional arguments affecting the fastCrr procedure.
#'
#' @details Calculates the CIF using \code{fcrr} output conditional on \code{cov}.
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
#' fit <- fastCrr(ftime, fstatus, cov)
#' cov2 <- rnorm(5)
#' getCIF(fit, cov2)
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.

predict.fcrr <- function(fit, cov, getBootstrapVariance = TRUE, B = 100,
                   type = "none",
                   alpha = 0.05, seed = 1991,
                   tL = NULL, tU = NULL, ...){

  ## Error checking
  if(class(fit) != "fcrr") {
    stop("Object 'fit' must be of class fcrr")
  }

  if(is.null(fit$breslowJump)) {
    stop("Breslow jumps were not calculated. Please re-run model with 'getBreslowJumps = TRUE'")
  }

  if(is.null(fit$df)) {
    stop("Ordered data frame not returned. Please re-run model with 'returnDataFrame = TRUE'")
  }

  if(!(type %in% c("none", "bands", "interval"))) {
    type = "none"
    warning("type is incorrectly specified. Valid options are 'bands', 'interval', 'none'.
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
    tU <- max(fit$uftime)
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

  res  <- data.frame(ftime = fit$uftime, CIF = CIF.hat, lower = NA, upper = NA)
  #Get SD via bootstrap
  set.seed(seed)
  if(getBootstrapVariance) {
    CIF.boot <- matrix(NA, nrow = B, ncol = length(fit$uftime))
    colnames(CIF.boot) <- round(fit$uftime, 3)
    ftime <- fit$df$ftime
    fstatus <- fit$df$fstatus
    n <- length(ftime)
    X <- as.matrix(fit$df[, -(1:2)])
    for(i in 1:B) {
      bsamp  <- sample(n, n, replace = TRUE) #Bootstrap sample index
      fit.bs <- fastCrr(ftime[bsamp], fstatus[bsamp], X[bsamp, ], variance = FALSE, ...)
      CIF.bs <- 1 - exp(-cumsum(exp(sum(cov * fit.bs$coef)) * fit.bs$breslowJump[, 2]))
      CIF.boot[i, ] <- evalstep(fit.bs$breslowJump$time,
                                stepf = CIF.bs,
                                subst = 1E-16,
                                newtime = fit$uftime)
    }
    rm(CIF.bs)
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
    } else if (type == "interval") {
      #If interval type if pointwise
      llim   <- CIF.hat + qnorm(1 - alpha / 2) * CIF.sd
      ulim   <- CIF.hat - qnorm(1 - alpha / 2) * CIF.sd
      res  <- data.frame(ftime = fit$uftime, CIF = exp(-exp(CIF.hat)), lower = exp(-exp(llim)), upper = exp(-exp(ulim)))
    }
  }
  #Subset corresponding to tL and tU
  res <- subset(res, res$ftime >= tL & res$ftime <= tU)
  class(res) <- "predict.fcrr"
  res$type <- type
  return(res)
}