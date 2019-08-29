#' Create a Competing Risk Object
#'
#' @description Create a competing risk object, used as a response variable in the model formula for \code{fastCrr} and \code{fastCrrp}.
#' Adapted from the \code{Surv} object.
#' @param ftime A vector of event/censoring times.
#' @param fstatus A vector with unique code for each event type and a separate code for censored observations.
#' @param failcode Integer: code of \code{fstatus} that event type of interest (default is 1)
#' @param cencode Integer: code of \code{fstatus} that denotes censored observations (default is 0)
#' @param silent Logical: print information about coding.
#'
#' @return Returns an object, used as a response variable, of class \code{Crisk}.
#' \item{time}{vector of observed event times}
#' \item{status}{vector of event indicators. 0 = censored, 1 = event of interest, 2 = competing risks}
#'
#' @importFrom survival Surv
#' @export
#' @examples
#' library(fastcmprsk)
#'
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' obj <- Crisk(ftime, fstatus, silent = FALSE)
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#'
#' @seealso Surv
Crisk <- function(ftime, fstatus, cencode = 0, failcode = 1, silent = TRUE) {

  # Check for errors
  if(!cencode %in% unique(fstatus)) stop("cencode must be a valid value from fstatus")
  if(!failcode %in% unique(fstatus)) stop("cencode must be a valid value from fstatus")
  if(any(ftime < 0)) stop("all values of ftime must be positive valued")

  crisk.ind <- setdiff(fstatus, c(cencode, failcode))
  if(!silent) {
    writeLines(paste("Observations with the following fstatus values are considered as right censored:", cencode))
    writeLines(paste("Observations with the following fstatus values are considered as observing the primary event:", failcode))

    if(identical(crisk.ind, integer(0))) {
      warning("Competing risks are not present")
    } else {
      writeLines(paste(c("Observations with the following fstatus values are considered as competing risks:", sort(crisk.ind)), collapse = " "))
    }
  }

  obj <- suppressWarnings(Surv(ftime, fstatus)) # Suppress warning given by Surv function
  obj[, 2] <- fstatus
  obj[obj[, 2] != cencode & obj[, 2] != failcode, 2] = 2 #Convert competing risks to 2
  class(obj) <- "Crisk"
  return(obj)
}
