#' Plots solution path for penalized methods
#'
#' @description  Plots solution path for penalized methods
#'
#' @param x \code{fcrrp} object (output from \code{fastCrrp()})
#' @param ... additional arguments to \code{plot()}
#' @details Plots solution path for penalized methods. x-axis: log tuning parameter values. y-axis: coeffcient estimates.
#' @import graphics
#' @export
#'

plot.fcrrp <-
  function(object, ...) {
    plot(NA, xlab = expression(log[10](lambda[n])), ylab = expression(beta[j]),
         ylim = c(min(object$coef), max(object$coef)),
         xlim = c(log10(min(object$lambda.path)), log10(max(object$lambda.path))),
         ...)
    for(i in 1:dim(object$coef)[1]) {
      lines(x$coef[i, ] ~ log10(object$lambda.path), col = i)
    }
  }
