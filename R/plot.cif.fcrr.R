#' Plots predicted cumulative incidence function
#'
#' @description  Plots predicted cumulative incidence function
#'
#' @param x \code{predict.fcrr} object (output from \code{predict(fcrr x)})
#' @param ... additional arguments to \code{plot()}
#' @import graphics
#' @return A plot of the estimated cumulative incidence function.
#' @export
#'
plot.predict.fcrr <-
  function(x, ...) {
    if(x$type == "none") {
      plot(x$CIF ~ x$ftime, xlab = "Time", ylab = "Estimated CIF",
           type = "s", ...)
    } else {
      plot(x$CIF ~ x$ftime, xlab = "Time", ylab = "Estimated CIF",
           ylim = c(min(x$lower), max(x$upper)),
           xlim = c(min(x$ftime), max(x$ftime)), type = "s", ...)
      lines(x$lower ~ x$ftime, lty = 2, type = "s")
      lines(x$upper ~ x$ftime, lty = 2, type = "s")
    }
  }
