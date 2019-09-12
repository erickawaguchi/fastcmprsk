#' Confidence Intervals for Model Parameters
#'
#' @description  Computes confidence intervals for one or more parameters in a fitted model of class \code{fcrr}.
#'
#' @param object \code{fcrr} object (output from \code{fastCrr()})
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required
#' @param digits Number of significant difits to round to.
#' @param ... Additional arguments. Not implemented.
#' @details Prints out table of confidence intervals for the Fine-Gray model.
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each coefficient estimate.
#' @export

confint.fcrr <-
  function(object, parm, level = 0.95, digits = max(options()$digits - 5, 2), ...) {

    if(!object$isVariance) {
      stop("Variance must be calculated. Rerun fastCrr with 'variance = TRUE'.")
    } else {
      ses <- sqrt(diag(vcov(object)))
    }

    cf <- coef(object)
    if(is.null(names(cf))) {
      pnames = paste0("x", 1:length(cf))
    } else {
      pnames = names(cf)
    }
    if (is.matrix(cf))
      cf <- setNames(as.vector(cf), pnames)
    if (missing(parm))
      parm <- pnames
    else if (is.numeric(parm))
      parm <- pnames[parm]
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    pct = paste0(a * 100, "%")
    ci = matrix(NA, nrow = length(cf), ncol = 2)
    ci[, 1] = cf + fac[1] * ses
    ci[, 2] = cf + fac[2] * ses
    rownames(ci) = pnames
    colnames(ci) = pct
    ci = ci[parm, ]
    ci
  }
