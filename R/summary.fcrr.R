#' Prints summary of a fcrr object
#'
#' @description  Prints summary statistics of a fcrr object
#'
#' @param x \code{fcrr} object (output from \code{fastCrr()})
#' @param ... additional arguments to \code{print()}
#' @details Prints the convergence status,
#' log-pseudo likelihood, the estimated coefficients, the estimated standard errors, and the two-sided p-values for the test of the individual coefficients equal to 0.
#' @export

summary.fcrr <-
  function(x,...) {
    cat('convergence:', x$converged, '\n')
    cat('# of iterations:', x$iter, '\n')

    cat('coefficients:\n')
    print(signif(x$coef,4), ...)
    if(!x$getVariance) {
      v <- x$var
    } else {
      v <- sqrt(diag(x$var))
    }
    cat('standard errors:\n')
    print(signif(v, 4),...)
    v <- 2 * (1 - pnorm(abs(x$coef) / v))
    cat('two-sided p-values:\n')
    print(signif(v, 2), ...)
    cat('log-pseudo likelihood:', x$logLik, '\n')
    invisible()
  }
