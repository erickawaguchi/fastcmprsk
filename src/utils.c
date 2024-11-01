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

// Get Euclidean norm
double getNorm(double *x, int p) {
  double x_norm = 0;
  for (int j = 0; j < p; j++) x_norm = x_norm + pow(x[j], 2);
  x_norm = sqrt(x_norm);
  return(x_norm);
}

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
  const int n = nin;
  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set
  double loglik = 0; //store loglik

  //Pointers to be freed later
  double *accSum = calloc(n, sizeof(double)); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) accSum[i] = 0;

  //Note: t2, ici, and x should be ordered in DECREASING order. (First time is largest)
  //Backward Scan [O(n)]
  for (int i = 0; i < n; i++)
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
  for(int i2 = (n - 1); i2 >= 0; i2--) {
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


  //calculate loglik [O(n)]
  for (int i = 0; i < n; i++) {
    if (ici[i] != 1) continue;
    loglik  -= log(accSum[i]);
  }
  free(accSum);
  return loglik;
}


//Calculate Breslow Jumps for Cumulative Hazard.
void getBreslowJumps(double *t2, int *ici, double *x, int *ncov, int *nin, double *wt, double *b, double *bj)
{
  // Look at Eq (2) from Fu et al. 2017.
  const int p = ncov[0],  n = nin[0];
  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set

  //Pointers to be freed later
  double *eta = calloc(n, sizeof(double)); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *accSum = calloc(n, sizeof(double)); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) accSum[i] = 0;

  for (int i = 0; i < n; i++) {
    //initialize values to 0
    eta[i] = 0;

    for (int j = 0; j < p; j++)
      eta[i] += b[j] * x[n * j + i];
  }

  //Note: t2, ici, and x should be ordered in DECREASING order. (First time is largest)
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
  for(int i2 = (n - 1); i2 >= 0; i2--) {
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


  int count = 0; // count number of events. Breslow jumps only occur at these event times
  //calculate loglik [O(n)]
  for (int i = 0; i < n; i++) {
    if (ici[i] != 1) continue;
    bj[count] = (1 / accSum[i]);
    count += 1;
  }

  free(eta);
  free(accSum);
}


//Calculate gradient and hessian diagonals for fixed beta
void getGradientAndHessian(double *t2, int *ici, int *nin, double *wt,
                           double *eta, double *st, double *w, double *lik)
{
  const int  n = nin[0];

  // pointers (to be freed later)
  double *accNum1 = calloc(n, sizeof(double)); //accumulate the backwards numerator
  for (int i = 0; i < n; i++) accNum1[i] = 0;
  double *accNum2 = calloc(n, sizeof(double)); //acumulate the foreward numerator (weighted)
  for (int i = 0; i < n; i++) accNum2[i] = 0;
  double *accSum = calloc(n, sizeof(double)); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) accSum[i] = 0;

  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set
  double loglik = 0;

  //initialization
  for (int i = 0; i < n; i++)
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
  for(int i2 = (n - 1); i2 >= 0; i2--) {
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
  tmp1 = 0; tmp2 = 0; //reset temporary vals

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
  tmp1 = 0;
  tmp2 = 0; //reset tmp vals
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

  for (int i = 0; i < n; i++) {
    if (ici[i] != 1) continue;
    loglik  -= log(accSum[i]);
  }

  *lik = loglik;
  free(accNum1);
  free(accNum2);
  free(accSum);
}

//Calculate null gradient
void getNullGradient(double *t2, int *ici, int *nin, double *wt, double *st)
{
  const int  n = nin[0];

  // pointers (to be freed later)
  double *accNum1 = calloc(n, sizeof(double)); //accumulate the backwards numerator
  for (int i = 0; i < n; i++) accNum1[i] = 0;
  double *accNum2 = calloc(n, sizeof(double)); //acumulate the foreward numerator (weighted)
  for (int i = 0; i < n; i++) accNum2[i] = 0;
  double *accSum = calloc(n, sizeof(double)); //accumulate sum over both accNum1 and accNum2
  for (int i = 0; i < n; i++) accSum[i] = 0;

  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set
  //initialization
  for (int i = 0; i < n; i++)
  {
    tmp1 += 1; //track cumulative sum over 1:n
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
      tmp2 += 1 / wt[i2];
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
  tmp1 = 0; tmp2 = 0; //reset temporary vals

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
    st[i] = accNum1[i];
  }

  //Perform linear scan for competing risks
  tmp1 = 0;
  tmp2 = 0; //reset tmp vals
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
    st[i] += accNum1[i] / wt[i];
  }

  for (int i = 0; i < n; i++) {
    if(ici[i] != 1) {
      st[i] = - st[i];
    } else {
      st[i] = (1 - st[i]);
    }
  }

  free(accNum1);
  free(accNum2);
  free(accSum);
}
