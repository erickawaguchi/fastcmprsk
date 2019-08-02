#' Extract variance-covariance matrix from an "fcrr" object.
#'
#' @description  Similar functional utility to \code{vcov} methods.
#'
#' @param object \code{fcrr} object.
#' @param ... Additional arguments. Not implemented.
#' @export
#'
vcov.fcrr <- function(object, ...) {
  if(!object$isVariance) {
    stop("Variance must be calculated. Rerun fastCrr with 'variance = TRUE'.")
  } else {
    object$var
  }
}
