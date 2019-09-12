#' Akaike's An Information Criterion
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrr} object
#' @param k Numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ... Additional arguments. Not implemented.
#' @return A numeric value with the corresponding AIC (or BIC, or ..., depending on \code{k}).
#' @export
#'
AIC.fcrr <- function(object, ..., k = 2)
  -2 * object$logLik + k * sum(coef(object) != 0)
