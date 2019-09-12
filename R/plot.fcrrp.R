#' Plots solution path for penalized methods
#'
#' @description  Plots solution path for penalized methods
#'
#' @param x \code{fcrrp} object (output from \code{fastCrrp()})
#' @param ... additional arguments to \code{plot()}
#' @details Plots solution path for penalized methods. x-axis: log tuning parameter values. y-axis: coeffcient estimates.
#' @import graphics
#' @return A plot of the solution path for the chosen penalized method.
#' @export
#'

plot.fcrrp <-
  function(x, ...) {
    plot(NA, xlab = expression(log[10](lambda[n])), ylab = expression(beta[j]),
         ylim = c(min(x$coef), max(x$coef)),
         xlim = c(log10(min(x$lambda.path)), log10(max(x$lambda.path))),
         ...)
    for(i in 1:dim(x$coef)[1]) {
      lines(x$coef[i, ] ~ log10(x$lambda.path), col = i)
    }
  }
