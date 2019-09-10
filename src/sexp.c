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
// SEXP arguments that are used across several functions


//////////////////////////////////////////////////////////////////////////////////////
//Standardize design matrix
// Returns standardized design matrix, vector of column means, vector of column sd.
SEXP standardize(SEXP X_) {
    // Declarations
    int n = nrows(X_);
    int p = ncols(X_);
    SEXP XX_, c_, s_;
    PROTECT(XX_ = allocMatrix(REALSXP, n, p));
    PROTECT(c_ = allocVector(REALSXP, p));
    PROTECT(s_ = allocVector(REALSXP, p));
    double *X = REAL(X_);
    double *XX = REAL(XX_);
    double *c = REAL(c_);
    double *s = REAL(s_);

    for (int j = 0; j < p; j++) {

        // Center (Calculate mean and subtract)
        c[j] = 0;
        for (int i = 0; i < n; i++) {
            c[j] += X[j * n + i];
        }
        c[j] = c[j] / n;
        for (int i = 0; i < n; i++) XX[j * n + i] = X[j * n + i] - c[j];

        // Scale (Calculate sdev and divide)
        s[j] = 0;
        for (int i = 0; i < n; i++) {
            s[j] += pow(XX[j * n + i], 2);
        }
        s[j] = sqrt(s[j] / n);
        for (int i = 0; i < n; i++) XX[j * n + i] = XX[j * n + i] / s[j];
    }

    // Return list
    SEXP res;
    PROTECT(res = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(res, 0, XX_); // Standardized design matrix
    SET_VECTOR_ELT(res, 1, c_); // Mean
    SET_VECTOR_ELT(res, 2, s_); // Standard deviations
    UNPROTECT(4);
    return(res);
}

//////////////////////////////////////////////////////////////////////////////////////
// Wrapper for evaluating log-likelihood for different beta
SEXP evalLogLikelihood(SEXP t2_, SEXP ici_, SEXP eta_, SEXP wt_) {
    //Declaration
    int n = length(t2_);
    //int L = length(lambda);

    // initialize
    double *eta = REAL(eta_);
    double *t2 = REAL(t2_);
    double *wt = REAL(wt_);
    int *ici = INTEGER(ici_);
    double loglik = getLogLikelihood(t2, ici, eta, wt, n);
    return(ScalarReal(loglik));
}

//////////////////////////////////////////////////////////////////////////////////////
SEXP getResultsCrr(SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian, SEXP linpred) {
    SEXP res;
    PROTECT(res = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(res, 0, beta); //coefficient estimates
    SET_VECTOR_ELT(res, 1, Dev); //deviance = -2*loglik
    SET_VECTOR_ELT(res, 2, iter); //iterations until convergence
    SET_VECTOR_ELT(res, 3, residuals); //residuals
    SET_VECTOR_ELT(res, 4, score); //gradient
    SET_VECTOR_ELT(res, 5, hessian); //hessian
    SET_VECTOR_ELT(res, 6, linpred); //hessian
    UNPROTECT(8);
    return(res);
}

SEXP getResultsCrrp(SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian, SEXP linpred, SEXP converged) {
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 8));
  SET_VECTOR_ELT(res, 0, beta); //coefficient estimates
  SET_VECTOR_ELT(res, 1, Dev); //deviance = -2*loglik
  SET_VECTOR_ELT(res, 2, iter); //iterations until convergence
  SET_VECTOR_ELT(res, 3, residuals); //residuals
  SET_VECTOR_ELT(res, 4, score); //gradient
  SET_VECTOR_ELT(res, 5, hessian); //hessian
  SET_VECTOR_ELT(res, 6, linpred); //hessian
  SET_VECTOR_ELT(res, 7, converged); //check convergence
  UNPROTECT(9);
  return(res);
}

