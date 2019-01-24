#ifndef SPARSE_UTILS_H
#define SPARSE_UTILS_H

// This is the content of the .h file, which is where the declarations go
void getCoordinateMatrixVector(double *x, int *row_idx, int *col_idx,
                               int *ncov, double *y, double *val);
double getSparseWeightedCrossProduct(double *x, int *row_ptr, int *col_idx, double *y, double *w, int n);
double getSparseWeightedSumSquares(double *x, int *row_ptr, int *col_idx, double *w, int n);
// This is the end of the header guard
#endif
