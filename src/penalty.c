#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#include "utils.h"
#define LEN sizeof(double)

//////////////////////////////////////////////////////////////////////////////////////
// ERIC S. KAWAGUCHI
//////////////////////////////////////////////////////////////////////////////////////
// Utilities for Package
// Create penalty function for ridge and lasso regressions.


double getRidgeDelta(double grad, double hess, double a, double lam) {
    return (-(grad + a * lam) / (hess + lam));
}

double getLassoDelta(double grad, double hess, double a, double lam) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= lam) return(-a);
  else return(s * (fabs(z) - lam) / (hess) - a);
}

double getMcpDelta(double grad, double hess, double a, double lam, double gamma) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= lam) return(-a);
  else if (fabs(z) <= gamma * lam) return(s * (fabs(z) - lam)/(hess * (1 - 1 / gamma)) - a);
  else return(z / hess - a);
}

double getScadDelta(double grad, double hess, double a, double lam, double gamma) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= lam) return(-a);
  else if (fabs(z) <= 2 * lam) return(s * (fabs(z) - lam) / hess - a);
  else if (fabs(z) <= gamma * lam) return(s * (fabs(z) - gamma * lam / (gamma - 1)) / (hess * (1 - 1 / (gamma - 1))) - a);
  else return(z / hess - a);
}
