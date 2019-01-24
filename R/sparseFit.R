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

sparseCrr <- function(outcomes, covariates, failcode = 1, cencode = 0,
                    eps = 1E-6,
                    max.iter = 1000, getBreslowJumps = TRUE,
                    startingCoefficients = NULL,
                    getVariance = TRUE,
                    varianceControl = varianceControl(),
                    quiet = FALSE){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")
  if(!is.data.frame(outcomes)) stop("outcomes must be a data frame.")
  if(!is.data.frame(covariates)) stop("covariates must be a data frame.")
  #outcomes = dat$outcomes[, c(1, 8, 2)]

  if(!quiet) writeLines("Preprocessing data prior to model fit.")

  #Reorder observations in covariate data frame
  n <- dim(outcomes)[1]
  outcomes <- outcomes[order(outcomes[, 2], -outcomes[, 3], decreasing = TRUE), ]
  covariates[, 1] = rep(order(outcomes[ , 1]), tabulate(covariates[, 1]))

  #Sort by covariate
  covariates = covariates[order(covariates[, 2]), ]
  row_ptr = c(0, cumsum(tabulate(covariates[, 2]))) #Where next covariates start in the COO format

  cencode = 0; failcode = 1
  cenind         <- ifelse(outcomes[, 3] == cencode, 1, 0)
  outcomes[, 3]  <- ifelse(outcomes[, 3] == failcode, 1, 2 * (1 - cenind))
  u <- do.call('survfit', list(formula = Surv(outcomes[, 2], cenind) ~ 1))

  # uuu is weight function (IPCW)
  u <- approx(c(0, u$time, max(u$time) * (1 + 10 * .Machine$double.eps)), c(1, u$surv, 0),
              xout = outcomes[, 2] * (1 - 100 * .Machine$double.eps), method = 'constant',
              f = 0, rule = 2)
  uuu <- u$y

  p = length(row_ptr) - 1

  if(!quiet) writeLines("Fitting model.")
  sparseFit   <- .Call("ccd_sparse", covariates[, 3],
                       as.integer(row_ptr), as.integer(covariates[, 1] - 1),
                       as.integer(p),
                       as.numeric(outcomes[, 2]), as.integer(outcomes[, 3]), uuu,
                       eps, as.integer(max.iter), PACKAGE = "fastcmprsk")[1:3]

  if (sparseFit[[3]] == max.iter) {
    warning("Maximum number of iterations reached. Estimates may not be stable")
  }

  bhat = sparseFit[[1]]

  val <- structure(list(coef = bhat,
                        logLik = sparseFit[[2]] / -2,
                        iter = sparseFit[[3]],
                        breslowJump = getBreslowJumps,
                        uftime = unique(rev(outcomes[, 2][outcomes[, 3] == 1])),
                        call = sys.call()),
                   class = "fcrr")

  val
}
