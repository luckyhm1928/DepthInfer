#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Call calls */
  extern SEXP _DepthInfer_AtoB(SEXP);
extern SEXP _DepthInfer_fofrMRtest(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DepthInfer_fofrMRtestBTS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DepthInfer_FourierBasis(SEXP, SEXP);
extern SEXP _DepthInfer_FPCAcov(SEXP, SEXP, SEXP);
extern SEXP _DepthInfer_kFoldCVcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _DepthInfer_orthoL2Equidense(SEXP, SEXP);
extern SEXP _DepthInfer_rowMins(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_DepthInfer_AtoB",             (DL_FUNC) &_DepthInfer_AtoB,             1},
  {"_DepthInfer_fofrMRtest",       (DL_FUNC) &_DepthInfer_fofrMRtest,       8},
  {"_DepthInfer_fofrMRtestBTS",    (DL_FUNC) &_DepthInfer_fofrMRtestBTS,    9},
  {"_DepthInfer_FourierBasis",     (DL_FUNC) &_DepthInfer_FourierBasis,     2},
  {"_DepthInfer_FPCAcov",          (DL_FUNC) &_DepthInfer_FPCAcov,          3},
  {"_DepthInfer_kFoldCVcpp",       (DL_FUNC) &_DepthInfer_kFoldCVcpp,       4},
  {"_DepthInfer_orthoL2Equidense", (DL_FUNC) &_DepthInfer_orthoL2Equidense, 2},
  {"_DepthInfer_rowMins",          (DL_FUNC) &_DepthInfer_rowMins,          1},
  {NULL, NULL, 0}
};

void R_init_DepthInfer(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}