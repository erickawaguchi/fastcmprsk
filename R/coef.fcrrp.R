#' Extract coefficients from an "fcrrp" object.
#'
#' @description  Similar functional utility to \code{coef} methods.
#'
#' @param object \code{fcrrp} object
#' @param ... Additional arguments. Not implemented.
#' @return Coefficients extracted from the model object \code{object}.
#' @export
#'
coef.fcrrp <- function(object, ...)
  object$coef
