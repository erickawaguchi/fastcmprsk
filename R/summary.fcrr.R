#' Summary method for fastCrr
#'
#' @description  Generate and print summaries of \code{fastCrr} output.
#'
#' @param object \code{fcrr} x (output from \code{fastCrr()})
#' @param conf.int Logical. Whether or not to outut confidence intervals.
#' @param alpha Significance level of the confidence intervals.
#' @param digits Numer of significant difits to round to.
#' @param ... additional arguments to \code{print()}
#' @details The summary method produces an ANOVA table for the coefficient estimates of the Fine-Gray model.
#' @return The form of the value returned by \code{summary} depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.
#' @export

summary.fcrr <-
  function(object, conf.int = TRUE, alpha = 0.05, digits = max(options()$digits - 5, 2), ...) {

    if(!object$isVariance) {
      se <- NA
    } else {
      se <- sqrt(diag(object$var))
    }

    beta <- object$coef
    out <- list(call = object$call, converged = object$converged,
                iterations = object$iter,
                logLik = object$logLik,
                logLik.null = object$logLik.null,
                ncov = length(object$coef),
                lrt = object$lrt)
    if(is.null(names(beta))) {
      names(beta) <- paste0("x", 1:length(beta))
    }
    tmp <- cbind(beta, exp(beta), se, beta / se,
                 signif(2 * (1 - pnorm(abs(beta) / se)), digits))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                         "se(coef)", "z value", "Pr(>|z|)"))

    out$coef <- tmp
    if(conf.int)
    { a <- alpha / 2
    a <- c(a, 1 - a)
    z <- qnorm(a)
    tmp <- cbind(exp(beta), exp(-beta),
                 exp(beta + z[1] * se), exp(beta + z[2] * se))
    dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                         paste(format(100 * a, trim = TRUE,
                                                      scientific = FALSE,
                                                      digits = 3), "%", sep="")))
    out$conf.int <- tmp
    }
    class(out) <- "summary.fcrr"
    out
  }
