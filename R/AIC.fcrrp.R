#' Akaike's An Information Criterion
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrrp} object
#' @param k Numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ... Additional arguments. Not implemented.
#' @export
#'
AIC.fcrrp <- function(object, ..., k = 2)
  -2 * object$logLik + k * apply(object$coef, 2, function(x) sum(x != 0))
