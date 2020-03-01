#' Controls for Variance Calculation
#'
#' @description Controls for variance calculation for the fastcmprsk package.
#'
#' @param B Integer: Number of bootstrap samples needed for variance estimation.
#' @param seed Integer: Seed value for bootstrapping. Results may differ if parallelized.
#' @param useMultipleCores Logical: Set to TRUE if parallelizing. (Default is FALSE).
#' @param extractMatrix Logical: Extract matrix of bootstrap estimates (Default is FALSE)
#' @return Returns a list for variance options inputted into \code{fastCrr}.
#' \item{B}{same as what is defined in argument.}
#' \item{seed}{same as what is defined in argument.}
#' \item{mcores}{same as what is defined in argument \code{useMultipleCores}.}
#' \item{extract}{same as what is defined in argument \code{extractMatrix}.}
#'
#' @export
#' @details Variance-covariance estimation is done via bootstrap.
#' Independent bootstrap runs can be performed both in serial and parallel. Parallelization is done via the
#' \code{doParallel} package.
#' @examples
#'
#' library(fastcmprsk)
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' vc <- varianceControl(B = 100, seed = 2019, useMultipleCores = FALSE)
#' fit1 <- fastCrr(Crisk(ftime, fstatus) ~ cov, variance = TRUE, var.control = vc)
#' fit1$var # Estimated covariance matrix via bootstrap
#'

varianceControl <- function(B = 100L, seed = 1991L, useMultipleCores = FALSE, extractMatrix = FALSE)
{

  if (B <= 0) {
    warning("The value of 'B' must be a non-negative integer. Set to 100")
    B <- 100L
  }

  if(seed <= 0) {
    warnings("The value of 'seed' must be a non-negative integer. Set to 1991")
    seed <- 1991L
  }


  obj          <- list()
  obj$B        <- B
  obj$seed     <- seed
  obj$mcores   <- useMultipleCores
  obj$extract  <- extractMatrix
  return(obj)
}
