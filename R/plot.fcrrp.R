#' Plots solution path for penalized methods
#'
#' @description  Plots solution path for penalized methods
#'
#' @param x \code{fcrrp} object (output from \code{fastCrrp()})
#' @param ... additional arguments to \code{plot()}
#' @details
#' @export
#'

plot.fcrrp <-
  function(x, ...) {
    plot(NA, xlab = "log(Lambda)", ylab = "beta",
         ylim = c(min(x$coef), max(x$coef)),
         xlim = c(log(min(x$lambda.path)), log(max(x$lambda.path))),
         ...)
    for(i in 1:dim(x$coef)[1]) {
      lines(x$coef[i, ] ~ log(x$lambda.path), col = i)
    }
  }
