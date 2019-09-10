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
#define LEN sizeof(double)


//////////////////////////////////////////////////////////////////////////////////////
//Coordinate Descent Algorithm for Dense Fine-Gray Model
SEXP ccd_dense(SEXP x_, SEXP t2_, SEXP ici_, SEXP wt_,
               SEXP eps_, SEXP max_iter_) {

    //Declaration
    int n = length(t2_);
    int p = length(x_) / n;

    //Output
    SEXP res, beta, Dev, iter, residuals, score, hessian, linpred;
    PROTECT(beta = allocVector(REALSXP, p));
    double *b = REAL(beta);
    for (int j = 0; j < p; j++) b[j] = 0;
    PROTECT (score = allocVector(REALSXP, n));
    double *s = REAL(score);
    for (int i = 0; i < n; i++) s[i] = 0;
    PROTECT (hessian = allocVector(REALSXP, n));
    double *h = REAL(hessian);
    for (int i = 0; i < n; i++) h[i] = 0;
    PROTECT(residuals = allocVector(REALSXP, n));
    double *r = REAL(residuals);
    PROTECT(Dev = allocVector(REALSXP, 2));
    for (int i = 0; i < 2; i++) REAL(Dev)[i] = 0;
    PROTECT(iter = allocVector(INTSXP, 1));
    for (int i = 0; i < 1; i++) INTEGER(iter)[i] = 0;
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

    //internal storage
    double nullDev; //to store null deviance
    double grad, hess, shift,  si, delta;
    int converged; //convergence check
    double tmp1 = 0; //track backward sum for uncensored events risk set
    double tmp2 = 0; //track forward sum for competing risks risk set
    //end of declaration;

    //Start regression
    //calculate null deviance
    nullDev = -2 * getLogLikelihood(t2, ici, eta, wt, n);
    REAL(Dev)[0] = nullDev; // save null deviance

    //start
    while (INTEGER(iter)[0] < max_iter) {
        if (REAL(Dev)[1] - nullDev > 0.99 * nullDev) break;

        INTEGER(iter)[0]++;

        //Reset values
        tmp1 = 0, tmp2 = 0;
        for (int i = 0; i < n; i++) {
            st[i] = 0;
            w[i] = 0;
        }

        //Backward Scan [O(n)]
        for (int i = 0; i < n; i++)
        {
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

        // calculate xwr and xwx & update beta_j
        for (int j = 0; j < p; j++) {
            grad = -getWeightedCrossProduct(x, w, r, n, j); // jth component of gradient [l'(b)]
            hess = getWeightedSumSquares(x, w, n, j); // jth component of hessian [l''(b)]
            delta = -grad / hess;

            // Do one dimensional ridge update.
            // Employ trust region as in Genkin et al. (2007) for quadratic approximation.
            b[j] = a[j] + sgn(delta) * fmin(fabs(delta), diffBeta[j]);
            diffBeta[j] = fmax(2 * fabs(delta), diffBeta[j] / 2);

            // Update r
            shift = b[j] - a[j];
            if (shift != 0) {
                for (int i = 0; i < n; i++) {
                    si = shift * x[j * n + i];
                    r[i] -= si;
                    eta[i] += si;
                }
            } //end shift
        } //for j = 0 to (p - 1)

        // Check for convergence
        converged = checkConvergence(b, a, eps, 0, p);
        for (int j = 0; j < p; j++)
            a[j] = b[j];

        //Calculate deviance
        REAL(Dev)[1] = -2 * getLogLikelihood(t2, ici, eta, wt, n);

        for (int i = 0; i < n; i++){
            lp[i] = eta[i];
            s[i] = st[i];
            h[i] = w[i];
        }
        if (converged)  break;
    } //for while loop

    // Free Calloc variables:
    Free(a);
    Free(eta);
    Free(st);
    Free(w);
    Free(diffBeta);
    Free(accNum1);
    Free(accNum2);
    Free(accSum);

    res = getResultsCrr(beta, Dev, iter, residuals, score, hessian, linpred);
    return(res);
}


//////////////////////////////////////////////////////////////////////////////////////

