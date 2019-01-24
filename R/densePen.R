#' Penalized Fine-Gray Model Estimation
#'
#' @description Performs penalized regression for the proportional subdistribution hazards model.
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
#' @param getBootstrapVariance Logical: Get standard error estimates for parameter estimates via bootstrap
#'
#' @details The \code{pshBAR} function penalizes the log-partial likelihood of the proportional subdistribution hazards model
#' from Fine and Gray (1999) with the Broken Adaptive Ridge (BAR) penalty. A cyclic coordinate descent algorithm is used for implementation.
#' For stability, the covariate matrix \code{X} is standardized prior to implementation.
#'
#' Special cases: Fixing \code{xi} and \code{lambda} to 0 results in the standard competing risk regression using \code{crr}.
#' Fixing \code{lambda} to 0 and specifying \code{xi} will result in a ridge regression solution.
#' @return Returns a list of class \code{pshBAR}.
#'
#' @import survival doParallel
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

fastCrrp <- function(ftime, fstatus, X, failcode = 1, cencode = 0,
                    eps = 1E-6,
                    max.iter = 1000, getBreslowJumps = TRUE,
                    penalty = c("lasso", "ridge", "mcp", "scad"),
                    lambda = NULL,
                    penalty.factor = rep(1, ncol(X)),
                    gamma = switch(penalty, scad = 3.7, 2.7)){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")
  if(!(penalty %in% c("lasso", "ridge", "mcp", "scad"))) stop("penalty is incorrectly specified. Please select lasso, ridge, mcp, or scad")
  if(min(lambda) < 0) stop("lambda must be a non-negative number.")
  if (gamma <= 1 & penalty == "mcp")
    stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty == "scad")
    stop("gamma must be greater than 2 for the SCAD penalty")
  # Sort time
  n <- length(ftime)
  p <- ncol(X)
  d <- data.frame(ftime = ftime, fstatus = fstatus)
  if (!missing(X)) d$X <- as.matrix(X)
  d        <- d[order(d$ftime, -d$fstatus, decreasing = TRUE), ]
  ftime    <- d$ftime
  cenind   <- ifelse(d$fstatus == cencode, 1, 0)
  fstatus  <- ifelse(d$fstatus == failcode, 1, 2 * (1 - cenind))
  X <- d$X
  u <- do.call('survfit', list(formula = Surv(ftime, cenind) ~ 1,
                               data = data.frame(ftime, cenind)))

  # uuu is weight function (IPCW)
  u <- approx(c(0, u$time, max(u$time) * (1 + 10 * .Machine$double.eps)), c(1, u$surv, 0),
              xout = ftime * (1 - 100 * .Machine$double.eps), method = 'constant',
              f = 0, rule = 2)
  uuu <- u$y

  # Standardize design matrix here
  std    <- .Call("standardize", X, PACKAGE = "fastcmprsk")
  XX     <- std[[1]]
  center <- std[[2]]
  scale  <- std[[3]]
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ , nz, drop = FALSE]


  # Order lambda in decreasing order increasing order. [Dense -> Sparse Model]
  lambda <- sort(lambda)

  # Fit the PSH penalized model
  denseFit   <- .Call("ccd_dense_pen", XX, as.numeric(ftime), as.integer(fstatus), uuu,
                      eps, as.integer(max.iter), penalty, as.double(lambda),
                      as.double(penalty.factor), as.double(gamma), PACKAGE = "fastcmprsk")

  bhat <- matrix(denseFit[[1]], p, length(lambda)) / scale
  colnames(bhat) <- round(lambda, 4)
  # Calculate Breslow Baseline
  if(getBreslowJumps) {
    jump = matrix(NA, ncol = length(lambda) + 1, nrow = length(unique(ftime[fstatus == 1])))
    jump[, 1] = time = unique(rev(ftime[fstatus == 1]))
    for(l in 1:length(lambda)) {
    bjump = .C("getBreslowJumps", as.double(ftime), as.integer(fstatus), as.double(X),
               as.integer(p), as.integer(n), as.double(uuu), as.double(bhat[, l]), double(sum(fstatus == 1)),
               PACKAGE = "fastcmprsk")
    jump[, l + 1] = as.vector(rev(unique(bjump[[8]])) * table(ftime[fstatus == 1], fstatus[fstatus == 1]))
    }
    colnames(jump) = c("time", paste0("Lam:", round(lambda, 4)))
    getBreslowJumps <- data.frame(jump)
  } #End Breslow jump

  #Do not calculate bootstrap variances for penalized model. Can be unstable


  #Results to store:
  val <- structure(list(coef = bhat,
                        logLik = denseFit[[2]][-1] / -2,
                        logLik.null = denseFit[[2]][1] / -2,
                        iter = denseFit[[3]],
                        converged = denseFit[[8]],
                        breslowJump = getBreslowJumps,
                        uftime = unique(rev(ftime[fstatus == 1])),
                        gamma = gamma,
                        call = sys.call()),
                   class = "fcrr")
  val
}
