#' Extract coefficients from an "fcrr" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrr} object
#' @param ... Additional arguments. Not implemented.
#' @return Coefficients extracted from the model object \code{object}.
#' @export
#'
coef.fcrr <- function(object, ...)
  object$coef
