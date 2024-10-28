#' Penalized Fine-Gray Model Estimation via two-way linear scan
#'
#' @description Performs penalized regression for the proportional subdistribution hazards model.
#' Penalties currently include LASSO, MCP, SCAD, and ridge regression. User-specificed weights can be assigned
#' to the penalty for each coefficient (e.g. implementing adaptive LASSO and broken adaptive ridge regerssion).
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a Crisk object as returned by the \code{Crisk} function.
#' @param data a data.frame in which to interpret the variables named in the formula.
#' @param eps Numeric: algorithm stops when the relative change in any coefficient is less than \code{eps} (default is \code{1E-6})
#' @param max.iter Numeric: maximum iterations to achieve convergence (default is 1000)
#' @param getBreslowJumps Logical: Output jumps in Breslow estimator for the cumulative hazard (one for each value of lambda).
#' @param standardize Logical: Standardize design matrix.
#' @param penalty Character: Penalty to be applied to the model. Options are "lasso", "scad", "ridge", "mcp", and "enet".
#' @param lambda A user-specified sequence of \code{lambda} values for tuning parameters.
#' @param alpha L1/L2 weight for elastic net regression.
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of \code{lambda.max} (if \code{lambda} is NULL).
#' @param nlambda Number of \code{lambda} values (default is 25).
#' @param penalty.factor A vector of weights applied to the penalty for each coefficient. Vector must be of length equal to the number of columns in \code{X}.
#' @param gamma Tuning parameter for the MCP/SCAD penalty. Default is 2.7 for MCP and 3.7 for SCAD and should be left unchanged.
#'
#' @details The \code{fastCrrp} functions performed penalized Fine-Gray regression.
#' Parameter estimation is performed via cyclic coordinate descent and using a two-way linear scan approach to efficiently
#' calculate the gradient and Hessian values. Current implementation includes LASSO, SCAD, MCP, and ridge regression.
#' @return Returns a list of class \code{fcrrp}.
#' \item{coef}{fitted coefficients matrix with \code{nlambda} columns and \code{nvars} columns}
#' \item{logLik}{vector of log-pseudo likelihood at the estimated regression coefficients}
#' \item{logLik.null}{log-pseudo likelihood when the regression coefficients are 0}
#' \item{lambda.path}{sequence of tuning parameter values}
#' \item{iter}{number of iterations needed until convergence at each tuning parameter value}
#' \item{converged}{convergence status at each tuning parameter value}
#' \item{breslowJump}{Jumps in the Breslow baseline cumulative hazard (used by \code{predict.fcrr})}
#' \item{uftime}{vector of unique failure (event) times}
#' \item{penalty}{same as above}
#' \item{gamma}{same as above}
#' \item{above}{same as above}
#'
#' @import survival
#' @export
#' @examples
#'
#' library(fastcmprsk)
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit <- fastCrrp(Crisk(ftime, fstatus) ~ cov, lambda = 1, penalty = "RIDGE")
#' fit$coef
#'
#' @references
#' Fu, Z., Parikh, C.R., Zhou, B. (2017) Penalized variable selection in competing risks
#' regression. \emph{Lifetime Data Analysis} 23:353-376.
#'
#' Breheny, P. and Huang, J. (2011) Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection. \emph{Ann. Appl. Statist.}, 5: 232-253.
#'
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#'
#' Kawaguchi, E.S., Shen J.I., Suchard, M. A., Li, G. (2020) Scalable Algorithms for Large Competing Risks Data, Journal of Computational and Graphical Statistics

fastCrrp <- function(formula, data,
                    eps = 1E-6,
                    max.iter = 1000,
                    getBreslowJumps = TRUE,
                    standardize = TRUE,
                    penalty = c("LASSO", "RIDGE", "MCP", "SCAD", "ENET"),
                    lambda = NULL, alpha = 0,
                    lambda.min.ratio = 0.001, nlambda = 25,
                    penalty.factor,
                    gamma = switch(penalty, scad = 3.7, 2.7)){

  ## Error checking
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")
  if(!(penalty %in% c("LASSO", "RIDGE", "MCP", "SCAD", "ENET")))
    stop("penalty is incorrectly specified. Please select 'LASSO', 'RIDGE', 'MCP', 'SCAD' or 'ENET'.")
  if (gamma <= 1 & penalty == "MCP")
    stop("gamma must be greater than 1 for the MCP penalty")
  if (gamma <= 2 & penalty == "SCAD")
    stop("gamma must be greater than 2 for the SCAD penalty")
  if(alpha < 0 | alpha > 1) stop("alpha must be between 0 and 1")

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

  # Create data-driven lambda path if one is not provided
  if(is.null(lambda)) {
    if(lambda.min.ratio < 0 | lambda.min.ratio > 1) stop("lambda.min.ratio must be between 0 and 1.")
    if(nlambda < 1) stop("nlambda must be larger than one")

    eta0 <- rep(0, n) #Linear predictor when beta = 0
    sw <- .C("getGradientAndHessian", as.double(dat$ftime), as.integer(dat$fstatus),
             as.integer(n), as.double(dat$wt),
             as.double(eta0), double(n), double(n), double(1),
             PACKAGE = "fastcmprsk") #Linearized version of crrp function
    score0 <- sw[[6]]
    w0 <- sw[[8]]
    r0 <- ifelse(w0 == 0, 0, score0 / w0)
    z <- eta0 + r0
    lambda.max <- max(t(w0 * z) %*% dat$X) / n # TO DO: Import this into C
    lambda = 10^(seq(log10(lambda.max), log10(lambda.min.ratio * lambda.max), len = nlambda))
  }



  # Order lambda in decreasing order increasing order. [Dense -> Sparse Model]
  if(min(lambda) < 0) stop("lambda(s) must be non negative.")
  lambda  <- sort(lambda, decreasing = TRUE)
  nlambda <- length(lambda)

  # Fit the PSH penalized model
  if(missing(penalty.factor)) penalty.factor = rep(1, p)

  if(penalty %in% c("LASSO", "RIDGE", "MCP", "SCAD")) {
  denseFit   <- .Call("ccd_dense_pen", dat$X, as.numeric(dat$ftime), as.integer(dat$fstatus), dat$wt,
                      eps, as.integer(max.iter), penalty, as.double(lambda),
                      as.double(penalty.factor), as.double(gamma), PACKAGE = "fastcmprsk")

  bhat <- matrix(denseFit[[1]], p, length(lambda)) / dat$scale
  } else {
    #Elastic Net regression
    denseFit   <- .Call("ccd_dense_enet", dat$X, as.numeric(dat$ftime), as.integer(dat$fstatus), dat$wt,
                        eps, as.integer(max.iter), as.double(alpha), as.double(lambda),
                        as.double(penalty.factor), PACKAGE = "fastcmprsk")
    bhat <- matrix(denseFit[[1]], p, length(lambda)) / dat$scale
  }
  colnames(bhat) <- round(lambda, 4)

  # Calculate Breslow Baseline
  if(getBreslowJumps) {
    jump = matrix(NA, ncol = length(lambda) + 1, nrow = length(unique(ftime[fstatus == 1])))
    jump[, 1] = unique(rev(dat$ftime[dat$fstatus == 1]))
    for(l in 1:length(lambda)) {
      bjump = .C("getBreslowJumps", as.double(dat$ftime), as.integer(dat$fstatus), as.double(sweep(sweep(dat$X, 2, dat$scale, "*"), 2, dat$center, `+`)),
                 as.integer(p), as.integer(n), as.double(dat$wt), as.double(bhat[, l]), double(sum(fstatus == 1)),
                 PACKAGE = "fastcmprsk")
      jump[, l + 1] = as.vector(rev(unique(bjump[[8]])) * table(dat$ftime[dat$fstatus == 1], dat$fstatus[dat$fstatus == 1]))
    }
    colnames(jump) = c("time", paste0("Lam = ", round(lambda, 4)))
    getBreslowJumps <- data.frame(jump)
  } #End Breslow jump

  #Results to store:
  val <- structure(list(coef = bhat,
                        logLik = denseFit[[2]][-1] / -2,
                        logLik.null = denseFit[[2]][1] / -2,
                        lambda.path = lambda,
                        iter = denseFit[[3]],
                        converged = denseFit[[8]],
                        breslowJump = getBreslowJumps,
                        uftime = unique(rev(dat$ftime[dat$fstatus == 1])),
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        call = sys.call()),
                   class = "fcrrp")
  val
}
