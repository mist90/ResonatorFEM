#include "MathEighValues.h"
#include <vector>
#include <cstddef>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <Eigen/Sparse>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

/* Generalized symmetric-definite eigenproblem  A x = lambda B x  solved densely
 * with LAPACK's dsygv.
 *
 * The original code used ARPACK++/SuperLU shift-invert, but the vendored
 * ARPACK++ 1.2 headers call an obsolete SuperLU ABI (dgstrf lost `drop_tol` and
 * gained a GlobalLU_t* in SuperLU >=5), which corrupts the stack at runtime.
 * For the cavity resonator the free-edge system is modest in size, symmetric,
 * and B (the epsilon-weighted mass matrix) is positive definite, so a dense
 * LAPACK solve is both robust and exact — it returns the *entire* spectrum,
 * including the large edge-element null space near lambda=0.
 *
 * Returns eigenvalues in ascending order and the corresponding eigenvectors
 * packed column-major: eigVector[row + mode*n] (same layout MicroEngine expects
 * for rootsGlobalMatrix). For very large meshes a sparse shift-invert solver
 * (Spectra/SLEPc) should replace this — see docs/PLAN.md. */

extern "C" {
void dsygv_(const int* itype, const char* jobz, const char* uplo,
            const int* n, double* a, const int* lda, double* b, const int* ldb,
            double* w, double* work, const int* lwork, int* info);
}

bool MathEighValVector(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB,
                       std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.width() != matrixA.height()) return false;
    if(matrixB.width() != matrixB.height()) return false;
    if(matrixA.width() != matrixB.width())  return false;

    const int n = (int)matrixA.width();
    if(n <= 0) return false;

    /* Dense column-major copies (A x = lambda B x). Both matrices are symmetric
     * and fully populated by the assembly, so we read every element. */
    std::vector<double> A((std::size_t)n * n, 0.0);
    std::vector<double> B((std::size_t)n * n, 0.0);
    for(int col = 0; col < n; col++)
        for(int row = 0; row < n; row++)
        {
            A[(std::size_t)row + (std::size_t)col * n] = matrixA.element((uint32_t)col, (uint32_t)row);
            B[(std::size_t)row + (std::size_t)col * n] = matrixB.element((uint32_t)col, (uint32_t)row);
        }

    eighValue.assign(n, 0.0);
    const int itype = 1;          /* A x = lambda B x            */
    const char jobz = 'V';        /* eigenvalues and eigenvectors */
    const char uplo = 'U';
    int lda = n, ldb = n, info = 0, lwork = -1;

    /* Workspace query. */
    double workQuery = 0.0;
    dsygv_(&itype, &jobz, &uplo, &n, A.data(), &lda, B.data(), &ldb,
           eighValue.data(), &workQuery, &lwork, &info);
    if(info != 0) return false;

    lwork = (int)workQuery;
    if(lwork < 1) lwork = 1;
    std::vector<double> work((std::size_t)lwork, 0.0);

    dsygv_(&itype, &jobz, &uplo, &n, A.data(), &lda, B.data(), &ldb,
           eighValue.data(), work.data(), &lwork, &info);
    if(info != 0) {
        std::fprintf(stderr, "dsygv failed: n=%d info=%d (%s)\n", n, info,
                     info > n ? "B (mass) not positive definite" : "no convergence");
        return false;   /* info>n: B not positive definite; else no convergence */
    }

    /* On success A holds the eigenvectors as columns (column-major):
     * eigVector[row + mode*n]. */
    eigVector = A;
    return true;
}

/* Convert a symmetric, fully-populated MathMatrixSparse to Eigen sparse. */
static Eigen::SparseMatrix<double> toEigen(MathMatrixSparse<double>& M, int n)
{
    std::vector<Eigen::Triplet<double> > trip;
    for(int col = 0; col < n; col++)
        for(int row = 0; row < n; row++)
        {
            double v = M.element((uint32_t)col, (uint32_t)row);
            if(v != 0.0) trip.emplace_back(row, col, v);
        }
    Eigen::SparseMatrix<double> E(n, n);
    E.setFromTriplets(trip.begin(), trip.end());
    return E;
}

/* Core shift-invert: `nev` eigenpairs of Ae x = lambda Be x nearest `sigma`.
 * Retries perturbed shifts if (Ae - shift*Be) hits an eigenvalue. */
static bool shiftInvertCore(const Eigen::SparseMatrix<double>& Ae,
                            const Eigen::SparseMatrix<double>& Be,
                            double sigma, int nev,
                            Eigen::VectorXd& vals, Eigen::MatrixXd& vecs)
{
    const int n = (int)Ae.rows();
    if(nev < 1) nev = 1;
    if(nev > n - 2) nev = n - 2;
    int ncv = std::min(n, std::max(4 * nev + 1, 60));

    using OpType  = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    BOpType Bop(Be);

    const double perturb[] = { 1.0, 1.02, 0.98, 1.05, 0.95, 1.09, 0.91 };
    for(double f : perturb) {
        double s = sigma * f;
        try {
            OpType op(Ae, Be);
            Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
                geigs(op, Bop, nev, ncv, s);
            geigs.init();
            geigs.compute(Spectra::SortRule::LargestMagn, 2000, 1e-10);   /* nearest s */
            if(geigs.info() == Spectra::CompInfo::Successful) {
                vals = geigs.eigenvalues();
                vecs = geigs.eigenvectors();
                return true;
            }
        } catch (const std::exception&) {
            /* singular (Ae - s*Be); try the next perturbed shift. */
        }
    }
    std::fprintf(stderr, "shift-invert failed near sigma=%g (nev=%d n=%d)\n", sigma, nev, n);
    return false;
}

/* Pack Eigen results into the (eighValue, eigVector[row + mode*n]) layout. */
static void packResults(const Eigen::VectorXd& vals, const Eigen::MatrixXd& vecs, int n,
                        std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    const int k = (int)vals.size();
    eighValue.resize(k);
    for(int i = 0; i < k; i++) eighValue[i] = vals[i];
    eigVector.assign((std::size_t)k * n, 0.0);
    for(int mode = 0; mode < k; mode++)
        for(int row = 0; row < n; row++)
            eigVector[(std::size_t)row + (std::size_t)mode * n] = vecs(row, mode);
}

bool MathEighValVectorShiftInvert(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB,
                                  double sigma, int nev,
                                  std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.width() != matrixA.height()) return false;
    if(matrixB.width() != matrixB.height()) return false;
    if(matrixA.width() != matrixB.width())  return false;
    const int n = (int)matrixA.width();
    if(n <= 2) return false;

    Eigen::SparseMatrix<double> Ae = toEigen(matrixA, n);
    Eigen::SparseMatrix<double> Be = toEigen(matrixB, n);
    Eigen::VectorXd vals;
    Eigen::MatrixXd vecs;
    if(!shiftInvertCore(Ae, Be, sigma, nev, vals, vecs)) return false;
    packResults(vals, vecs, n, eighValue, eigVector);
    return true;
}

bool MathEighValVectorShiftInvertGauged(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB,
                                        const std::vector<std::pair<uint32_t, uint32_t> >& dofNodes,
                                        const std::vector<double>& dofLen,
                                        uint32_t nNodes, double penaltyS, double sigma, int nev,
                                        std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.width() != matrixA.height()) return false;
    if(matrixB.width() != matrixB.height()) return false;
    if(matrixA.width() != matrixB.width())  return false;
    const int n = (int)matrixA.width();
    if(n <= 2 || (int)dofNodes.size() != n || (int)dofLen.size() != n || nNodes == 0) return false;

    Eigen::SparseMatrix<double> Ae = toEigen(matrixA, n);   /* stiffness S */
    Eigen::SparseMatrix<double> Be = toEigen(matrixB, n);   /* mass T      */

    /* Discrete gradient G (n edges x nNodes): edge d = (tail a, head b), scaled
     * by edge length for this len*Whitney basis.
     *
     * EXPERIMENTAL / WIP: the exact discrete gradient of this hand-rolled,
     * length-scaled edge basis is not yet fully identified — empirically the
     * residual ||S*G||/||S|| bottoms out around 0.15 (not ~0), so G is not quite
     * in the null space of S and the penalty still perturbs the physical modes.
     * Consequently the grad-div penalty is OFF by default (penaltyFactor = 0);
     * getting the null-space projection exact is a focused follow-up. The
     * *formulation* below is correct given a correct G: P = (T G) D^-1 (Gᵀ T) is
     * zero on physical modes (Gᵀ T e = 0) and lifts the gradient modes. */
    std::vector<Eigen::Triplet<double> > gtrip;
    gtrip.reserve((std::size_t)n * 2);
    for(int d = 0; d < n; d++) {
        uint32_t a = dofNodes[d].first, b = dofNodes[d].second;
        if(a == UINT32_MAX || dofLen[d] <= 0.0) continue;
        double w = dofLen[d];
        gtrip.emplace_back(d, (int)a, -w);
        gtrip.emplace_back(d, (int)b, +w);
    }
    Eigen::SparseMatrix<double> G((int)n, (int)nNodes);
    G.setFromTriplets(gtrip.begin(), gtrip.end());

    /* Mass-metric grad-div penalty: P = (T G) D^-1 (G^T T), D = diag(G^T T G).
     * P e = 0 for physical modes (G^T T e = 0), so they are preserved exactly;
     * gradient (null-space) modes are lifted to ~penaltyS. */
    Eigen::SparseMatrix<double> TG = Be * G;                 /* n x nNodes    */
    Eigen::SparseMatrix<double> L  = Eigen::SparseMatrix<double>(G.transpose()) * TG; /* nNodes^2 = G^T T G */
    Eigen::VectorXd dinv(nNodes);
    for(uint32_t k = 0; k < nNodes; k++) {
        double dk = L.coeff((int)k, (int)k);
        dinv[k] = (dk > 1e-30) ? 1.0 / dk : 0.0;
    }
    Eigen::SparseMatrix<double> P = TG * dinv.asDiagonal() * Eigen::SparseMatrix<double>(TG.transpose());
    Eigen::SparseMatrix<double> Aeff = Ae + penaltyS * P;

    Eigen::VectorXd vals;
    Eigen::MatrixXd vecs;
    if(!shiftInvertCore(Aeff, Be, sigma, nev, vals, vecs)) return false;
    packResults(vals, vecs, n, eighValue, eigVector);
    return true;
}
