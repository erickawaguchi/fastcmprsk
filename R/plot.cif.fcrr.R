#' Plots predicted cumulative incidence function
#'
#' @description  Plots predicted cumulative incidence function
#'
#' @param x \code{predict.fcrr} object (output from \code{predict(fcrr object)})
#' @param ... additional arguments to \code{plot()}
#' @import graphics
#' @export
#'
plot.predict.fcrr <-
  function(object, ...) {
    if(object$type == "none") {
      plot(object$CIF ~ object$ftime, xlab = "Time", ylab = "Estimated CIF",
           type = "s", ...)
    } else {
      plot(object$CIF ~ object$ftime, xlab = "Time", ylab = "Estimated CIF",
           ylim = c(min(object$lower), max(object$upper)),
           xlim = c(min(object$ftime), max(object$ftime)), type = "s", ...)
      lines(object$lower ~ object$ftime, lty = 2, type = "s")
      lines(object$upper ~ object$ftime, lty = 2, type = "s")
    }
  }
