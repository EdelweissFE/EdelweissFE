#pragma once

namespace PardisoInterface{
    int solve(const double* values, const int* innerIndices, const int* outerIndexPtr, int nnz, int cols, int rows,  const double* rhs, int nRhs, double* out);
}
