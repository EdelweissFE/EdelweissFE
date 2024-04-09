#include "klu.h"

int solve(double *LHS, int *innerIndices, int *outerIndexPtr, int n,
          double *rhs, int nRhs) {

  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common Common;
  klu_defaults(&Common);

  Symbolic = klu_analyze(n, outerIndexPtr, innerIndices, &Common);
  Numeric = klu_factor(outerIndexPtr, innerIndices, LHS, Symbolic, &Common);

  klu_solve(Symbolic, Numeric, n, nRhs, rhs, &Common);
  return 0;
}
