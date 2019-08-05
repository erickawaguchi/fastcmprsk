#' Extract log-pseudo likelihood from an "fcrr" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrr} object
#' @param ... Additional arguments. Not implemented.
#' @export
#'
logLik.fcrr <- function(object, ...)
  object$logLik
