#' Extract log-pseudo likelihood from an "fcrrp" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrrp} object
#' @param ... Additional arguments. Not implemented.
#' @return Returns the log-pseudo likelihood of object \code{object}.
#' @export
#'
logLik.fcrrp <- function(object, ...)
  object$logLik
