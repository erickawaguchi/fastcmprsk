#' Extract coefficients from an "fcrr" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrr} object
#' @param ... some methods for this generic function require additional arguments. This package does not.
#' @export
#'
coef.fcrr <- function(object, ...)
  object$coef
