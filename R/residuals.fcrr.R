#' Calculate residuals for 'fcrr' fit
#'
#' @description Calculates martingale and Cox-Snell residuaks for the Fine-Gray model.
#'
#' @param object Output from \code{fcrr} object.
#' @param type Character string indicating the type of residual desired. Possible values are \code{"cox-snell"} or \code{"martingale"}
#' which are based under a similar derivation to the Cox proportional hazards model.
#' @param ... additional arguments (currently unused).
#'
#' @details Calculates residuals. \code{returnDataFrame = TRUE} must be specified in \code{fastCrr}.
#'
#' @import survival dynpred
#' @export
#' @useDynLib fastcmprsk
#' @examples
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' fit <- fastCrr(Crisk(ftime, fstatus) ~ cov, returnDataFrame = TRUE)
#' residuals(fit, type = "martingale")
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.

residuals.fcrr <- function(object, type = "martingale", ...){

  ## Error checking
  if(class(object) != "fcrr") {
    stop("Object 'fit' must be of class fcrr")
  }

  if(is.null(object$breslowJump)) {
    stop("Breslow jumps were not calculated. Please re-run model with 'getBreslowJumps = TRUE'")
  }

  if(is.null(object$df)) {
    stop("Ordered data frame not returned. Please re-run model with 'returnDataFrame = TRUE'")
  }

  if(!(type %in% c("cox-snell", "martingale"))) {
    type = "martingale"
    warning("type is incorrectly specified. Valid options are 'cox-snell', or 'martingale'.
            Set to 'martingale'")
  }

  # Calculate Breslow-type cumulative baseline hazard:
  H <- evalstep(object$breslowJump[, 1],
                stepf = cumsum(object$breslowJump[, 2]),
                subst = 0,
                newtime = object$df$ftime)

  # Calculate exponentiated linear predictor (currently in R, may be slow for large datasets)
  ncov <- dim(object$df)[2]
  eeta <- exp(as.matrix(object$df[, 3:ncov]) %*% object$coef)

  rj  <- eeta * H # H_0 * exp(Zbeta)
  out <- list()

  if(type == "cox-snell") {
    out$res = rj
  } else {
    res = (object$df$fstatus == 1) - rj
  }

  out$residuals <- res
  out$ftime     <- object$df$ftime
  out$type      <- type
  return(out)
  }
