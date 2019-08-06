#Sorts data apporpriately and calcualtes IPCW
setupData <- function(ftime, fstatus, X, cencode, failcode, standardize) {
  d <- data.frame(ftime = ftime, fstatus = fstatus)
  if (!missing(X)) d$X <- as.matrix(X)
  idx      <- order(d$ftime, -d$fstatus, decreasing = TRUE)
  d        <- d[idx, ]
  ftime    <- d$ftime
  cenind   <- ifelse(d$fstatus == cencode, 1, 0)
  fstatus  <- ifelse(d$fstatus == failcode, 1, 2 * (1 - cenind))
  X <- d$X
  u <- do.call('survfit', list(formula = Surv(ftime, cenind) ~ 1,
                               data = data.frame(ftime, cenind)))

  # uuu is weight function (IPCW)
  u <- approx(c(0, u$time, max(u$time) * (1 + 10 * .Machine$double.eps)), c(1, u$surv, 0),
              xout = ftime * (1 - 100 * .Machine$double.eps), method = 'constant',
              f = 0, rule = 2)
  uuu <- u$y

  # Standardize design matrix here
  if(standardize) {
    std    <- .Call("standardize", X, PACKAGE = "fastcmprsk")
    XX     <- std[[1]]
    center <- std[[2]]
    scale  <- std[[3]]
    nz <- which(scale > 1e-6)
    if (length(nz) != ncol(XX)) XX <- XX[ , nz, drop = FALSE]
  } else {
    XX <- X
    scale <- 1
    center <- 0
  }
  return(list(
    ftime = ftime,
    fstatus = fstatus,
    X = XX,
    scale = scale,
    center = center,
    wt = uuu,
    idx = idx
  ))
}
