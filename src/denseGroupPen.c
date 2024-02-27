#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#include "utils.h"
#include "sexp.h"
#include "penalty.h"
#define LEN sizeof(double)


//////////////////////////////////////////////////////////////////////////////////////
// Group Coordinate Descent Algorithm for Dense Fine-Gray Model

SEXP ccd_dense_gpen(SEXP x_, SEXP t2_, SEXP ici_, SEXP K1_, SEXP wt_,
                    SEXP eps_, SEXP max_iter_, SEXP penalty_, SEXP lambda_, SEXP mult_,
                    SEXP gamma_) {

  //Declaration
  int n = length(t2_);
  int p = length(x_) / n;
  int L = length(lambda_);
  int J = length(K1_) - 1;
  int *K1 = INTEGER(K1_);

  //Output
  SEXP res, beta, Dev, iter, residuals, score, hessian, converged, linpred;
  PROTECT(beta = allocVector(REALSXP, L * p));
  double *b = REAL(beta);
  for (int j = 0; j < (L * p); j++) b[j] = 0;
  PROTECT (score = allocVector(REALSXP, L * n));
  double *s = REAL(score);
  for (int i = 0; i < (L * n); i++) s[i] = 0;
  PROTECT (hessian = allocVector(REALSXP, L * n));
  double *h = REAL(hessian);
  for (int i = 0; i <  (L * n); i++) h[i] = 0;
  PROTECT(residuals = allocVector(REALSXP, n));
  double *r = REAL(residuals);
  PROTECT(Dev = allocVector(REALSXP, L + 1));
  for (int i = 0; i < (L + 1); i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i = 0; i < L; i++) INTEGER(iter)[i] = 0;
  PROTECT(converged = allocVector(INTSXP, L));
  for (int i = 0; i < L; i++) INTEGER(converged)[i] = 0;
  PROTECT(linpred = allocVector(REALSXP, n));
  double *lp = REAL(linpred);
  for (int i = 0; i <  n; i++) lp[i] = 0;

  //Intermediate quantities for internal use (must be freed afterwards!)
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0;
  double *st = Calloc(n, double);
  for (int i = 0; i < n; i++) st[i] = 0;
  double *w = Calloc(n, double);
  for ( int i = 0; i < n; i++) w[i] = 0;
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *diffBeta = Calloc(p, double);
  for (int j = 0; j < p; j++) diffBeta[j] = 1;
  double *accNum1 = Calloc(n, double); //accumulate the backwards numerator
  for (int i = 0; i < n; i++) accNum1[i] = 0;
  double *accNum2 = Calloc(n, double); //acumulate the foreward numerator (weighted)
  for (int i = 0; i < n; i++) accNum2[i] = 0;
  double *accSum = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) accSum[i] = 0;


  //Pointers for R values
  double *x = REAL(x_);
  double *t2 = REAL(t2_);
  double *wt = REAL(wt_);
  int *ici = INTEGER(ici_);
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda_);
  double *m = REAL(mult_);
  double gamma = REAL(gamma_)[0];

  //internal storage
  double nullDev; //to store null deviance
  double grad, hess, shift, si, l1;
  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set
  //end of declaration;

  //Start regression
  //calculate null deviance
  nullDev = -2 * getLogLikelihood(t2, ici, eta, wt, n);
  REAL(Dev)[0] = nullDev;

  for(int l = 0; l < L; l++) {


    if (l != 0) {
      for (int j = 0; j < p; j++) {
        a[j] = b[(l - 1) * p + j]; //warm start
        //diffBeta[j] = 1; //reset diffBeta
      }
    }

    //start
    while (INTEGER(iter)[l] < max_iter) {
      if (REAL(Dev)[l + 1] - nullDev > 0.99 * nullDev) break;

      INTEGER(iter)[l]++; //increment iteration counter

      //Reset values
      tmp1 = 0, tmp2 = 0;
      //Backward Scan [O(n)]
      for (int i = 0; i < n; i++)
      {
        accSum[i] = 0;
        st[i] = 0;
        w[i] = 0;
        tmp1 += exp(eta[i]); //track cumulative sum over 1:n
        if (ici[i] != 1) {
          // if subject is not an event then accNum[i] = 0;
          accSum[i] = 0;
        } else {
          accSum[i] = tmp1;
        }
      }

      //Forward Scan (To take into account the competing risks component) [O(n)]
      for (int i2 = (n - 1); i2 >= 0; i2--) {
        if (ici[i2] == 2) {
          tmp2 += exp(eta[i2]) / wt[i2];
        }
        if (ici[i2] != 1) continue;
        accSum[i2] += wt[i2] * tmp2;
      }


      //taking into account ties [O(n)]
      for (int i2 = (n - 1); i2 >= 0; i2--) {
        if(ici[i2] == 2 || i2 == 0 || ici[i2 - 1] != 1) continue;
        if(t2[i2] == t2[i2 - 1]) {
          accSum[i2 - 1] = accSum[i2];
        }
      }


      //calculate score and hessian here
      double tmp1 = 0; tmp2 = 0; //reset temporary vals

      //linear scan for non-competing risks (backwards scan)
      for (int i = (n - 1); i >= 0; i--) {
        if(ici[i] == 1) {
          tmp1 += 1 / accSum[i];
          tmp2 += 1 / pow(accSum[i], 2);
          accNum1[i] = tmp1;
          accNum2[i] = tmp2;
        } else {
          accNum1[i] = tmp1;
          accNum2[i] = tmp2;
        }
      }

      //Fix ties here:
      for (int i = 0; i < n; i++) {
        //only needs to be adjusted consective event times
        if(ici[i] != 1 || i == (n - 1) || ici[i + 1] != 1) continue;
        if(t2[i] == t2[i + 1]) {
          accNum1[i + 1] = accNum1[i];
          accNum2[i + 1] = accNum2[i];
        }
      }


      //Store into st and w so we can reuse accNum1 and accNum2
      for (int i = 0; i < n; i++) {
        st[i] = accNum1[i] * exp(eta[i]);
        w[i] = accNum2[i] * pow(exp(eta[i]), 2);
      }

      //Perform linear scan for competing risks
      tmp1 = 0; tmp2 = 0; //reset tmp vals
      for (int i = 0; i < n; i++) {
        accNum1[i] = 0;
        accNum2[i] = 0;
        if(ici[i] == 1) {
          tmp1 += wt[i] / accSum[i];
          tmp2 += pow(wt[i] / accSum[i], 2);
        }
        if(ici[i] != 2) continue;
        accNum1[i] = tmp1;
        accNum2[i] = tmp2;
      }

      //Now combine to produce score and hessian
      for (int i = 0; i < n; i++) {
        //First, update st and w and then get score and hessian
        st[i] += accNum1[i] * (exp(eta[i]) / wt[i]);
        w[i] += accNum2[i] * pow(exp(eta[i]) / wt[i], 2);
      }

      for (int i = 0; i < n; i++) {
        w[i] = (st[i] - w[i]);
        if(ici[i] != 1) {
          st[i] = - st[i];
        } else {
          st[i] = (1 - st[i]);
        }
      }

      //end calculation of score and hessian

      for (int i = 0; i < n; i++){
        if (w[i] == 0) r[i] = 0;
        else r[i] = st[i] / w[i];
      }


      // Group Coordinate Update
      for (int g = 0; g < J; g++) {
        l1 = lam[l] * m[g] *sqrt(K1[g + 1] - K1[g]);

        // Calculate z
        int K = K1[g + 1] - K1[g];
        double *z = Calloc(K, double);
        double *v = Calloc(K, double);

        for (int j = K1[g]; j < K1[g + 1]; j++){
          z[j - K1[g]] = getWeightedCrossProduct(x, r, w, n, j) / n + getWeightedSumSquares(x, w, n, j) / n * a[j];
          v[j - K1[g]] = getWeightedSumSquares(x, w, n, j) / n;
        }

        double z_norm = getNorm(z, K);

        // Update b
        double len = 0;
        for (int j = K1[g]; j < K1[g + 1]; j++) {
          if (strcmp(penalty, "LASSO") == 0) len = getGroupLasso(z_norm, v[j - K1[g]], l1);
          if (strcmp(penalty, "MCP") == 0) len = getGroupMcp(z_norm, v[j - K1[g]],  l1, gamma);
          if (strcmp(penalty, "SCAD") == 0) len = getGroupScad(z_norm, v[j - K1[g]], l1, gamma);
          if ((len != 0) | (a[K1[g]] != 0)) {
            // If necessary, update b and r
            b[l * p + j] = len * z[j - K1[g]] / z_norm;
            double shift = b[l * p + j] - a[j];
            for (int i = 0; i < n; i++) {
              double si = shift * x[j * n + i];
              r[i] -= si;
              eta[i] += si;
            }
          }
        }
        // Free Calloc vars
        Free(z);
        Free(v);
      } // End group updating

      // Check for convergence
      INTEGER(converged)[l] = checkConvergence(b, a, eps, l, p);
      for (int j = 0; j < p; j++)
        a[j] = b[l * p + j];

      //Calculate deviance
      REAL(Dev)[l + 1] = -2 * getLogLikelihood(t2, ici, eta, wt, n);

      for (int i = 0; i < n; i++){
        lp[i] = eta[i];
        s[l * n + i] = st[i];
        h[l * n + i] = w[i];
      }
      if (INTEGER(converged)[l])  break;
    } //for while loop
  }

  //Free Calloc variables:
  Free(a);
  Free(eta);
  Free(st);
  Free(w);
  Free(diffBeta);
  Free(accNum1);
  Free(accNum2);
  Free(accSum);

  res =  getResultsCrrp(beta, Dev, iter, residuals, score, hessian, linpred, converged);
  return(res);
}


//////////////////////////////////////////////////////////////////////////////////////

