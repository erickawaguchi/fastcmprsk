#' Prints summary of a fcrr object
#'
#' @description  Prints summary statistics of a fcrr object
#'
#' @param x \code{fcrr} object (output from \code{fastCrr()})
#' @param ... additional arguments to \code{print()}
#' @details Prints the convergence status,
#' log-pseudo likelihood, the estimated coefficients, the estimated standard errors, and the two-sided p-values for the test of the individual coefficients equal to 0.
#' @export

print.summary.fcrr <- function (object, digits = max(options()$digits - 4, 3), ...)
{
  cat("Fine-Gray Regression via fastcmprsk package. \n\n")
  if(object$converged)
  { cat("fastCrr converged in", object$iterations, "iterations.\n \n")
  } else {
    cat("fastCrr did not converge. Estimates may be unstable.\n \n")
  }

  if(!is.null(object$call))
  { cat("Call:\n")
    dput(object$call)
    cat("\n")
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  print(object$coef)
  cat("\n")
  if(!is.null(object$conf.int))
    print(object$conf.int)

  cat("Pseudo Log-likelihood =", object$logLik, "\n")
  cat("Null Pseudo Log-likelihood =", object$logLik.null, "\n")
  cat("Pseudo likelihood ratio test =", object$lrt, "on", object$ncov, "df. \n")
  invisible()
}
