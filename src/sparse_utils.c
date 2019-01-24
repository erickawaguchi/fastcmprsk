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
//Coordinate Form Sparse Data Matrix Vector Algebra

// Calculate val = X * y for CSR data (needed to calculate loglikelihood)
void getSparseMatrixVector(double *x, int *row_ptr, int *col_idx,
                               int *nrows, double *y, double *val) {
  const int n = nrows[0];
  int i, k;
  for(i = 0; i < n; i++) {
    for(k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
      val[i] += x[k] * y[col_idx[k]];
    }
  }
}

// getWeightedCrossProduct and getWeightedSumSquares for CSR data
// For fixed i calculate val = x_i^T w * y
// Get x'diag(w)y = val
double getSparseWeightedCrossProduct(double *x, int *row_ptr, int *col_idx, double *y, double *w, int n) {
  double val = 0;
  for (int i = row_ptr[n]; i < row_ptr[n + 1]; i++)
    val += x[i] * y[col_idx[i]] * w[col_idx[i]];
  return(val);
}

// Get x'diag(w)x = val
double getSparseWeightedSumSquares(double *x, int *row_ptr, int *col_idx, double *w, int n) {
  double val = 0;
  for (int i = row_ptr[n]; i < row_ptr[n + 1]; i++)
    val += pow(x[i], 2) * w[col_idx[i]];
  return(val);
}

// Sparse versions of getLogLikelihood and getBreslowJumps
double getSparseLogLikelihood(double *t2, int *ici, int *row_ptr, int *col_idx, double *x, double *wt, double *b, int nin)
{
  int i, i2, k;
  const int n = nin;
  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set
  double loglik = 0; //store loglik

  //Pointers to be freed later
  double *eta = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (i = 0; i < n; i++) eta[i] = 0;
  double *accSum = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (i = 0; i < n; i++) accSum[i] = 0;

  for(i = 0; i < n; i++) {
      for(k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
        eta[i] += x[k] * b[col_idx[k]];
      }
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
void getSparseBreslowJumps(double *t2, int *ici, int *row_idx, int *col_idx, double *x, double *wt, double *b, int nin,
                             double *bj)
{
  int i, i2, k;
  const int n = nin;
  double tmp1 = 0; //track backward sum for uncensored events risk set
  double tmp2 = 0; //track forward sum for competing risks risk set

  //Pointers to be freed later
  double *eta = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (i = 0; i < n; i++) eta[i] = 0;
  double *accSum = Calloc(n, double); //accumulate sum over both accNum1 and accNum2
  for (i = 0; i < n; i++) accSum[i] = 0;

  //Current implementation requires O(nz) calculations to calculate linear predictor (Should try to optimize this later)
  for(i = 0; i < n; i++) {
    for(k = row_idx[i]; k < row_idx[i + 1]; k++) {
      eta[i] += x[k] * b[col_idx[k]];
    }
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
