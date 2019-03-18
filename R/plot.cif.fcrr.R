#' Plots predicted cumulative incidence function
#'
#' @description  Plots predicted cumulative incidence function
#'
#' @param x \code{cif.fcrr} object (output from \code{getCIF()})
#' @param ... additional arguments to \code{print()}
#' @export
#'
plot.cif.fcrr <-
  function(x, ...) {
    plot(x$CIF ~ x$ftime, xlab = "Time", ylab = "Estimated CIF",
         ylim = c(min(x$lower), max(x$upper)),
         xlim = c(min(x$ftime), max(x$ftime)), type = "s", ...)
    lines(x$lower ~ x$ftime, lty = 2,type = "s")
    lines(x$upper ~ x$ftime, lty = 2,type = "s")
  }
