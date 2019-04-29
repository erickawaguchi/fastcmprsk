#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#define LEN sizeof(double)

//////////////////////////////////////////////////////////////////////////////////////
// ERIC S. KAWAGUCHI
//////////////////////////////////////////////////////////////////////////////////////
// Utilities for Package

// Get x'y = val
double getMatrixVectorMultiplcation(double *X, double*y, int n, int j) {
    int nn = n * j;
    double val = 0;
    for (int i = 0; i < n; i++) val += X[nn + i] * y[i];
    return(val);
}
// Get x'diag(w)y = val
double getWeightedCrossProduct(double *X, double *y, double *w, int n, int j) {
    int nn = n * j;
    double val = 0;
    for (int i = 0; i < n; i++) val += X[nn + i] * y[i] * w[i];
    return(val);
}

// Get x'diag(w)x = val
double getWeightedSumSquares(double *X, double *w, int n, int j) {
    int nn = n * j;
    double val = 0;
    for (int i = 0; i < n; i++) val += w[i] * pow(X[nn + i], 2);
    return(val);
}

//Define sgn function: sgn(z) = 1 if z > 0, -1 if z < 0, 0 if z = 0
double sgn(double z) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  return(s);
}

// Criterion for convergence: Max relative error should be less than eps
int checkConvergence(double *beta, double *beta_old, double eps, int l, int p) {
    int converged = 1;
    for (int j = 0; j < p; j++) {
        if (fabs((beta[l * p + j] - beta_old[j]) / beta_old[j]) > eps) {
            converged = 0;
            break;
        }
    }
    return(converged);
}

// Calculate Log-Partial Likelihood
double getLogLikelihood(double *t2, int *ici, double *eta, double *wt, int nin)
{
    int i, i2;
    const int n = nin;
    double tmp1 = 0; //track backward sum for uncensored events risk set
    double tmp2 = 0; //track forward sum for competing risks risk set
    double loglik = 0; //store loglik

    //Pointers to be freed later
    double *accSum = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
    for (int i = 0; i < n; i++) accSum[i] = 0;

    //Note: t2, ici, and x should be ordered in DECREASING order. (First time is largest)
    //Backward Scan [O(n)]
    for (i = 0; i < n; i++)
    {
        tmp1 += exp(eta[i]); //track cumulative sum over 1:n
        if (ici[i] != 1) {
            // if subject is not an event then accNum[i] = 0;
            accSum[i] = 0;
        } else {
            loglik += eta[i];
            accSum[i] = tmp1;
        }
    }

    //Forward Scan (To take into account the competing risks component) [O(n)]
    for(i2 = (n - 1); i2 >= 0; i2--) {
        if (ici[i2] == 2) {
            tmp2 += exp(eta[i2]) / wt[i2];
        }
        if (ici[i2] != 1) continue;
        accSum[i2] += wt[i2] * tmp2;
    }


    //taking into account ties [O(n)]
    for(i2 = (n - 1); i2 >= 0; i2--) {
        if(ici[i2] == 2 || ici[i2 - 1] != 1 || i2 == 0) continue;
        if(t2[i2] == t2[i2 - 1]) {
            accSum[i2 - 1] = accSum[i2];
        }
    }


    //calculate loglik [O(n)]
    for(i = 0; i < n; i++) {
        if (ici[i] != 1) continue;
        loglik  -= log(accSum[i]);
    }
    Free(accSum);
    return loglik;
}


//Calculate Breslow Jumps for Cumulative Hazard.
void getBreslowJumps(double *t2, int *ici, double *x, int *ncov, int *nin, double *wt, double *b, double *bj)
{
  // Look at Eq (2) from Fu et al. 2017.
  int i, j, i2;
  const int p = ncov[0],  n = nin[0];
  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set

  //Pointers to be freed later
  double *eta = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *accSum = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) accSum[i] = 0;

  for(i = 0; i < n; i++) {
    //initialize values to 0
    eta[i] = 0;

    for (j = 0; j < p; j++)
      eta[i] += b[j] * x[n * j + i];
  }

  //Note: t2, ici, and x should be ordered in DECREASING order. (First time is largest)
  //Backward Scan [O(n)]
  for (i = 0; i < n; i++)
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
  for(i2 = (n - 1); i2 >= 0; i2--) {
    if (ici[i2] == 2) {
      tmp2 += exp(eta[i2]) / wt[i2];
    }
    if (ici[i2] != 1) continue;
    accSum[i2] += wt[i2] * tmp2;
  }


  //taking into account ties [O(n)]
  for(i2 = (n - 1); i2 >= 0; i2--) {
    if(ici[i2] == 2 || ici[i2 - 1] != 1 || i2 == 0) continue;
    if(t2[i2] == t2[i2 - 1]) {
      accSum[i2 - 1] = accSum[i2];
    }
  }


  int count = 0; // count number of events. Breslow jumps only occur at these event times
  //calculate loglik [O(n)]
  for(i = 0; i < n; i++) {
    if (ici[i] != 1) continue;
    bj[count] = (1 / accSum[i]);
    count += 1;
  }

  free(eta);
  free(accSum);
}
