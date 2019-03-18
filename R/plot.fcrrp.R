#' Plots solution path for penalized methods
#'
#' @description  Plots solution path for penalized methods
#'
#' @param x \code{fcrrp} object (output from \code{fastCrrp()})
#' @details
#' @export

plot.fcrrp <-
  function(x) {
    plot(NA, main = paste0("Solution path for ", toupper(x$penalty), "-penalized regression"), xlab = "log(Lambda)", ylab = "beta",
         ylim = c(min(x$coef), max(x$coef)),
         xlim = c(log(min(x$lambda.path)), log(max(x$lambda.path))))
    for(i in 1:dim(x$coef)[1]) {
      lines(x$coef[i, ] ~ log(x$lambda.path), col = i)
    }
  }
