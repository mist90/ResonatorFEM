#ifndef MATHEIGHVALUES
#define MATHEIGHVALUES

#include <vector>
#include <utility>
#include <cstdint>
#include <array>
#include <Eigen/SparseCore>

/* Generalized symmetric eigenproblem  A x = lambda B x  (A, B assembled directly
 * as Eigen sparse matrices). Eigenvectors are returned column-major:
 * eigVector[row + mode*n]. */

/* Eigensolver selection for the resonator solve. */
enum EigenSolverKind
{
    EIGSOLVER_SHIFT_INVERT = 0,   /* Spectra shift-invert on the grad-div-penalized operator:
                                     direct factorization, robust, memory-heavy */
    EIGSOLVER_AME          = 1    /* hypre AME (LOBPCG + AMS preconditioner) on the raw
                                     curl-curl operator: iterative, low memory */
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
 * the gradient null-space to ~penaltyS so the physical modes become smallest. */
bool MathEighValVectorShiftInvertGauged(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                                        const std::vector<std::pair<uint32_t, uint32_t> >& dofNodes,
                                        const std::vector<double>& dofLen,
                                        const std::vector<char>& interiorNode,
                                        uint32_t nNodes, double penaltyS, double sigma, int nev,
                                        std::vector<double>& eighValue, std::vector<double>& eigVector);

/* hypre AME: the `nev` smallest nonzero (physical) eigenpairs of the curl-curl
 * pencil, preconditioned with AMS. A, B are the free-DOF stiffness/mass in the
 * length-scaled Whitney basis; dofNodes/dofLen give each free DOF's reference
 * orientation (node serials) and edge length, pecEdges the eliminated (PEC)
 * edges, nodeCoords the vertex coordinates indexed by node serial. The gradient
 * null space is projected out by AME itself, so no grad-div penalty is needed.
 * Needs an EigenSolverRuntime alive in the process. LOBPCG iteration limit and
 * tolerance can be overridden with env AME_ITERS / AME_TOL. */
bool MathEighValVectorAME(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                          const std::vector<std::pair<uint32_t, uint32_t> >& dofNodes,
                          const std::vector<double>& dofLen,
                          const std::vector<std::pair<uint32_t, uint32_t> >& pecEdges,
                          const std::vector<std::array<double, 3> >& nodeCoords,
                          int nev,
                          std::vector<double>& eighValue, std::vector<double>& eigVector);

/* Process-wide MPI + hypre initialization for the AME solver. Create exactly one
 * in main() before any solve (MPI must be finalized on the thread that
 * initialized it, so this cannot be done lazily from a worker thread). */
class EigenSolverRuntime
{
public:
    EigenSolverRuntime(int* argc, char*** argv);
    ~EigenSolverRuntime();
    EigenSolverRuntime(const EigenSolverRuntime&) = delete;
    EigenSolverRuntime& operator=(const EigenSolverRuntime&) = delete;
};

#endif // MATHEIGHVALUES
