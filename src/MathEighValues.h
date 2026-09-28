#ifndef MATHEIGHVALUES
#define MATHEIGHVALUES

#include <vector>
#include <utility>
#include <cstdint>
#include <Eigen/SparseCore>

/* Generalized symmetric eigenproblem  A x = lambda B x  (A, B assembled directly
 * as Eigen sparse matrices). Eigenvectors are returned column-major:
 * eigVector[row + mode*n]. */

/* Eigensolver used by the gauged (grad-div penalty) path. */
enum EigenSolverKind
{
    EIGSOLVER_SHIFT_INVERT = 0,   /* Spectra shift-invert: direct factorization, robust, memory-heavy */
    EIGSOLVER_LOBPCG       = 1    /* iterative low-memory LOBPCG; falls back to shift-invert on failure */
};

/* Dense LAPACK dsygv: full spectrum (small/well-conditioned problems). */
bool MathEighValVector(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                       std::vector<double>& eighValue, std::vector<double>& eigVector);

/* Sparse shift-invert (Spectra): the `nev` eigenvalues nearest `sigma`. */
bool MathEighValVectorShiftInvert(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                                  double sigma, int nev,
                                  std::vector<double>& eighValue, std::vector<double>& eigVector);

/* As above, plus a mass-metric grad-div penalty A += penaltyS * (T G) D^-1 (Gᵀ T),
 * where G is the discrete gradient given by dofNodes[dof] = (nodeSerial1,nodeSerial2)
 * (UINT32_MAX skipped) scaled by 1/dofLen, restricted to interior nodes
 * (interiorNode[serial] != 0), and D = diag(Gᵀ T G). Zero on physical modes; lifts
 * the gradient null-space to ~penaltyS so the physical modes become smallest.
 * `solver` picks the eigensolver (see EigenSolverKind). */
bool MathEighValVectorShiftInvertGauged(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                                        const std::vector<std::pair<uint32_t, uint32_t> >& dofNodes,
                                        const std::vector<double>& dofLen,
                                        const std::vector<char>& interiorNode,
                                        uint32_t nNodes, double penaltyS, double sigma, int nev,
                                        EigenSolverKind solver,
                                        std::vector<double>& eighValue, std::vector<double>& eigVector);

#endif // MATHEIGHVALUES
