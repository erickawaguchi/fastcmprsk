#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
  extern void getBreslowJumps(void *, void *, void *, void *, void *, void *, void *, void *);
extern void getGradientAndHessian(void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
  extern SEXP ccd_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ccd_dense_enet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ccd_dense_pen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP standardize(SEXP);

static const R_CMethodDef CEntries[] = {
  {"getBreslowJumps",       (DL_FUNC) &getBreslowJumps,       8},
  {"getGradientAndHessian", (DL_FUNC) &getGradientAndHessian, 8},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"ccd_dense",      (DL_FUNC) &ccd_dense,       6},
  {"ccd_dense_enet", (DL_FUNC) &ccd_dense_enet,  9},
  {"ccd_dense_pen",  (DL_FUNC) &ccd_dense_pen,  10},
  {"standardize",    (DL_FUNC) &standardize,     1},
  {NULL, NULL, 0}
};

void R_init_fastcmprsk(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
