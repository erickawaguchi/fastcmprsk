#' Controls for Variance Calculation
#'
#' @description Controls for variance calculation for the fastcmprsk package.
#'
#' @param B Integer: Number of bootstrap samples needed for variance estimation.
#' @param parallel Logical: Whether or not to parallelize independent bootstrap runs.
#' @param ncores Integer: Number of cores needed if \code{parallel = TRUE}. Default is one less than total number of cores.
#' @param seed Integer: Seed value for bootstrapping. Results may differ is \code{parallel = TRUE}.
#' @return Returns a list for variance options inputted into \code{fastCrr}.
#' @export
#' @examples
#' library(fastcmprsk)
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2, 200, replace = TRUE)
#' cov <- matrix(runif(1000), nrow = 200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
#' vc <- varianceControl(B = 100, parallel = FALSE, seed = 2019)
#' fit1 <- fastCrr(ftime, fstatus, cov, getVariance = TRUE, var.control = vc)

varianceControl <- function(B = 100L, parallel = FALSE, ncores = 1, seed = 1991L)
{

  if(parallel) {
    if (ncores < 0L || ncores > detectCores() || is.null(ncores)) {
      warning("The value of 'ncores' was either 0, larger than the number of cores available of NULL. Set to 1")
      ncores <- 1L
    }
  }

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
  obj$parallel <- parallel
  obj$ncores   <- ncores
  obj$seed     <- seed

  return(obj)
}
