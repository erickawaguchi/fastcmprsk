#ifndef SEXP_H
#define SEXP_H

// Store ALL SEXP in quantities here
// This is the content of the .h file, which is where the declarations go
SEXP standardize(SEXP X_);
SEXP evalLogLikelihood(SEXP t2_, SEXP ici_, SEXP eta_, SEXP wt_);
SEXP cleanupCRR(double *a, double *eta, double *st, double *w, double *diffBeta, double *accNum1, double *accNum2, double *accSum,
                SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian, SEXP linpred);
SEXP cleanupCRRP(double *a,  double *eta, double *st, double *w, double *diffBeta, double *accNum1, double *accNum2, double *accSum,
                 SEXP beta, SEXP Dev, SEXP iter, SEXP residuals, SEXP score, SEXP hessian, SEXP linpred, SEXP converged);

// This is the end of the header guard
#endif
