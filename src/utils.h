#ifndef UTILS_H
#define UTILS_H

// This is the content of the .h file, which is where the declarations go
double ridge(double z, double l1, double v);
double getMatrixVectorMultiplcation(double *X, double*y, int n, int j);
double getWeightedCrossProduct(double *X, double *y, double *w, int n, int j);
double getWeightedSumSquares(double *X, double *w, int n, int j);
double sgn(double z);
int checkConvergence(double *beta, double *beta_old, double eps, int l, int p);
double getLogLikelihood(double *t2, int *ici, double *eta, double *wt, int nin);
void getBreslowJumps(double *t2, int *ici, double *x, int *ncov, int *nin, double *wt, double *b, double *bj);
// This is the end of the header guard
#endif
