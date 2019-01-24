#' Simulate data from the Fine-Gray Model using a Sparse Design Matrix
#'
#' @description  Simulate data from the model proposed in Fine and Gray (1999) for two causes. Cause 1 is assumed
#' to be of primary importance. Design matrix is stored in coordinate format.
#'
#' @param nobs Integer: Number of observations in simulated dataset.
#' @param beta1 A vector of effect sizes for cause 1 of length ncovs
#' @param beta2 A vector of effect sizes for cause 2 of length ncovs
#' @param u.min Numeric: controls lower bound of censoring distribution where C ~ U(u.min, u.max)
#' @param u.max Numeric: controls upper bound of censoring distribution where C ~ U(u.min, u.max)
#' @param p Numeric: value between 0 and 1 which controls the mixture probability.
#' @param numCovsPerObs Numeric: value between 0 and 1 which controls the proportion of non-zero covariate values per observation.
#' @param parallel: Logical: Specifies whether to parallelize the construction of the design matrix.
#' @param ncores: Integer: Number of cores to use when \code{parallel} is TRUE.
#' @param seednum: Integer: For reproducibility.
#' @details The function simulates data according to the setup by Fine and Gray (1999). See their paper for more information.
#' The \code{outcomes} object stores the row (observation) ID, observed failure time, and censoring indicator, respectively.
#' The \code{covariates} object stores the row ID, column ID, and value for the non-zero entry in the design matrix, respectively. (e.g. If x_14 = 3, then rowId = 1, covariateId = 4, and covariateValue = 3)
#' @return Returns a list with \code{outcomes}, \code{covariates}.
#' @examples
#' set.seed(10)
#' nobs <- 500
#' beta1 <- c(0.5, 0.4, 0, 0, 0.35, 0, 0, 0.7)
#' beta2 <- -beta1
#' X <- matrix(rnorm(nobs * 8), nrow = nobs)
#' dat <- simulateTwoCauseFineGrayModel(nobs, beta1, beta2, X, u.min = 0, u.max = 1, p = 0.5)
#' @references
#' Fine J. and Gray R. (1999) A proportional hazards model for the subdistribution of a competing risk.  \emph{JASA} 94:496-509.
#' @import doParallel
#' @export

simulateSparseFineGrayModel <- function(nobs = 1000, beta1, beta2, p = 0.5,
                                   u.min = 0, u.max = 1, numCovsPerObs = 0.1,
                                   parallel = TRUE, ncores = NULL, seednum = 1234) {

  #nobs: Number of observation rows
  #ncovs: Number of covariates
  #trueCoef: True value of the non-zero simulated regression coefficients
  #numCovs: Expected number of covariates observed by each observation.
  #fixed: Non-zero covariate positions are fixed
  #censor.rate: Percentage of censoring
  ncovs <- length(beta1)
  effectSizes <- data.frame(covariateId = 1:ncovs, hr1 = exp(beta1), hr2 = exp(beta2))

  set.seed(seednum)
  numCovs <- numCovsPerObs * ncovs
  covPerObs <- rpois(nobs, numCovs)
  covPerObs[covPerObs > ncovs] <- ncovs #In case number of expected covarates is greater than the number of actual covariates.
  covPerObs     <- data.frame(covPerObs = covPerObs)
  totalCovCount <- sum(covPerObs$covPerObs)

  # ----- Parallel or not parallel
  writeLines(paste("Generating Sparse Matrix"))

  if(parallel == FALSE) {
    covData       <- data.frame(rowId = rep(0, totalCovCount), covariateId = rep(0, totalCovCount),
                                covariateValue = rep(1, totalCovCount))
    track <- 1
    for (i in 1:nrow(covPerObs)) {
      #For the ith covariate, sample the observations that observe the ith covariate
      n <- covPerObs$covPerObs[i] #Number of observations that observe the ith covariate
      if (n != 0){
        covData$covariateId[track:(track + n - 1)] <- sample.int(size = n, ncovs)
        covData$rowId[track:(track + n - 1)] <- i
        track = track + n
      }
      cat(paste0(round(i / nobs * 100, 3), "%", " completed", " \r "))
    }; rm(track)
  } else if(parallel == TRUE) {
    #- Detect cores and set number of cores. If null then use max - 1
    numCores <- detectCores()
    if(is.null(ncores)) {
      ncores <- numCores[1] - 1
    } else if(ncores > numCores) {
      stop(paste0("ncores must be less than or equal to ", numCores))
    } else {
      ncores <- ncores
    }

    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    #- Create sparse matrix in parallel:
    finalMatrix <- foreach(i = 1:nrow(covPerObs), .combine = rbind) %dopar% {
      n <- covPerObs$covPerObs[i]
      tempMatrix <- matrix(0, nrow = n, ncol = 3)
      tempMatrix[, 1] <- i
      tempMatrix[, 2] <- sample.int(size = n, ncovs)
      tempMatrix[, 3] <- 1
      tempMatrix
    }
    stopCluster(cl)
    covData <- as.data.frame(finalMatrix)
    names(covData) <- c("rowId", "covariateId", "covariateValue")
  }

  # -----
  # Generate which subjects fail from cause 1:
  writeLines(paste("Generating Cause Indicators and Survival Times"))
  outcomes       <- data.frame(rowId = 1:nobs, y = 0, cause = 0)
  obs_hr         <- aggregate(cbind(hr1, hr2) ~ rowId, data = merge(covData, effectSizes), prod) # Every subjects exp(Xbeta)
  outcomes$cause <- 1 + rbinom(nobs, 1, prob = (1 - p)^obs_hr$hr1)

  outcomes <- merge(outcomes, obs_hr, all.x = TRUE)
  outcomes$hr1[is.na(outcomes$hr1)] <- 1 #Give observations with no covariates a hr of 1.
  outcomes$hr2[is.na(outcomes$hr2)] <- 1 #Give observations with no covariates a hr of 1.

  # ---- Generating data from Cox Model
  writeLines(paste("Generating Survival Times"))
  outcomes$survTime <- NA

  #Simulate Cause 1 subjects from subdistribution and cause 2 from exponential
  nc1 <- which(outcomes$cause == 1)
  nc2 <- which(outcomes$cause == 2)
  u1  <- runif(length(nc1))

  outcomes$survTime[nc1] <- -log(1 - (1 - (1 - u1 * (1 - (1 - p)^outcomes$hr1[nc1]))^(1 / outcomes$hr1[nc1])) / p)
  outcomes$survTime[nc2] <- rexp(n = length(nc2), rate = outcomes$hr2[nc2])
  outcomes$censorTime    <- runif(nobs, min = u.min, max = u.max)


  outcomes$time    <- pmin(outcomes$survTime, outcomes$censorTime)
  outcomes$y       <- outcomes$survTime == outcomes$time
  outcomes$y       <- outcomes$y * outcomes$cause

  sparseness <- 1 - (totalCovCount / (nobs * ncovs)) #Calculates sparsity
  writeLines(paste("Sparseness =", sparseness * 100,"%"))

  # ---- Store List
  out <- list()
  out$outcomes   <- outcomes[, c("rowId", "time", "y")]
  out$covariate  <- covData
  out$sparseness <- sparseness
  out$parallel   <- parallel
  out$seed       <- seednum
  return(out)
}

