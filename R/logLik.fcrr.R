#' Extract log-pseudo likelihood from an "fcrr" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrr} object
#' @param ... Additional arguments. Not implemented.
#' @return Returns the log-pseudo likelihood of object \code{object}.
#' @export
#'
logLik.fcrr <- function(object, ...)
  object$logLik
