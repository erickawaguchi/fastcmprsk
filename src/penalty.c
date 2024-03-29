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


double getRidge(double grad, double hess, double a, double lam) {
    double z = hess * a - grad;
    return(z / (hess + lam));
}

double getLasso(double grad, double hess, double a, double lam) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= lam) return(0);
  else return(s * (fabs(z) - lam) / (hess));
}

double getMcp(double grad, double hess, double a, double lam, double gamma) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= lam) return(0);
  else if (fabs(z) <= gamma * lam) return(s * (fabs(z) - lam)/(hess * (1 - 1 / gamma)));
  else return(z / hess);
}

double getScad(double grad, double hess, double a, double lam, double gamma) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= lam) return(0);
  else if (fabs(z) <= 2 * lam) return(s * (fabs(z) - lam) / hess);
  else if (fabs(z) <= gamma * lam) return(s * (fabs(z) - gamma * lam / (gamma - 1)) / (hess * (1 - 1 / (gamma - 1))));
  else return(z / hess);
}

//Include Elastic Net (06/04/2019)
double getElasticNet(double grad, double hess, double a, double lam, double alpha) {
  double z = hess * a - grad;
  int s = sgn(z);
  if (fabs(z) <= (lam * alpha)) return(0);
  else return(s * (fabs(z) - (lam * alpha)) / (hess + lam * (1 - alpha)));
}

//Include Group penalties (04/23/2020)
double getGroupMcp(double z, double hess, double lam, double gamma) {
  int s = sgn(z);
  if (fabs(z) <= lam) return(0);
  else if (fabs(z) <= gamma * lam) return(s * (fabs(z) - lam) / (hess * (1 - 1 / gamma)));
  else return(z / hess);
}

double getGroupScad(double z, double hess, double lam, double gamma) {
  int s = sgn(z);
  if (fabs(z) <= lam) return(0);
  else if (fabs(z) <= 2 * lam) return(s * (fabs(z) - lam) / hess);
  else if (fabs(z) <= gamma * lam) return(s * (fabs(z) - gamma * lam / (gamma - 1)) / (hess * (1 - 1 / (gamma - 1))));
  else return(z / hess);
}

double getGroupLasso(double z, double hess, double lam) {
  int s = sgn(z);
  if (fabs(z) <= lam) return(0);
  else return(s * (fabs(z) - lam) / (hess));
}
