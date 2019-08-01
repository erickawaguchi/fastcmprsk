#' Extract log-pseudo likelihood from an "fcrr" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrr} object
#' @param ... some methods for this generic function require additional arguments. This package does not.
#' @export
#'
logLik.fcrr <- function(object, ...)
  object$logLik
