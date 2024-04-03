#include <Eigen/Dense>
#include <Eigen/PardisoSupport>

namespace PardisoInterface{
    int solve(const double* values, const int* innerIndices, const int* outerIndexPtr, int nnz, int cols, int rows,  const double* rhs, int nRhs, double* out)
    {
        using namespace Eigen;
        typedef SparseMatrix<double, RowMajor> SpMat;
        typedef PardisoLU   <SpMat > Solver;
        Map<const SpMat > A(rows,cols,nnz,outerIndexPtr, innerIndices,values);
    
        Map< MatrixXd >               x(out, rows, nRhs) ;
        Map< const MatrixXd >         b(rhs, rows, nRhs) ;
    
        Solver solver;
        
        solver.analyzePattern(A);
        solver.factorize(A);
    
         x = solver.solve(b);
    
        return 0;
    }
}
