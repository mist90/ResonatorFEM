#include "MathEighValues.h"
#include <vector>
#include <cstddef>
#include <cstdio>
#include <algorithm>
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

bool MathEighValVectorShiftInvert(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB,
                                  double sigma, int nev,
                                  std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.width() != matrixA.height()) return false;
    if(matrixB.width() != matrixB.height()) return false;
    if(matrixA.width() != matrixB.width())  return false;

    const int n = (int)matrixA.width();
    if(n <= 2) return false;

    /* Convert to Eigen sparse (both matrices are symmetric and fully populated). */
    std::vector<Eigen::Triplet<double> > tripA, tripB;
    for(int col = 0; col < n; col++)
        for(int row = 0; row < n; row++)
        {
            double a = matrixA.element((uint32_t)col, (uint32_t)row);
            double b = matrixB.element((uint32_t)col, (uint32_t)row);
            if(a != 0.0) tripA.emplace_back(row, col, a);
            if(b != 0.0) tripB.emplace_back(row, col, b);
        }
    Eigen::SparseMatrix<double> Ae(n, n), Be(n, n);
    Ae.setFromTriplets(tripA.begin(), tripA.end());
    Be.setFromTriplets(tripB.begin(), tripB.end());

    if(nev < 1) nev = 1;
    if(nev > n - 2) nev = n - 2;
    /* A generous Arnoldi subspace helps convergence to interior eigenvalues near
     * a large null space. */
    int ncv = std::min(n, std::max(4 * nev + 1, 60));

    using OpType  = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    BOpType Bop(Be);

    /* Factorizing (A - shift*B) fails if `shift` coincides with an eigenvalue
     * (including a spurious one). Retry with small multiplicative perturbations
     * around the requested sigma until one factorizes and converges. */
    const double perturb[] = { 1.0, 1.02, 0.98, 1.05, 0.95, 1.09, 0.91 };
    Eigen::VectorXd vals;
    Eigen::MatrixXd vecs;
    bool ok = false;
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
                ok = true;
                break;
            }
            std::fprintf(stderr, "  shift %g: not converged\n", s);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "  shift %g: %s\n", s, e.what());
        }
    }
    if(!ok) {
        std::fprintf(stderr, "shift-invert failed near sigma=%g (nev=%d ncv=%d n=%d)\n",
                     sigma, nev, ncv, n);
        return false;
    }

    const int k = (int)vals.size();
    eighValue.resize(k);
    for(int i = 0; i < k; i++) eighValue[i] = vals[i];
    eigVector.assign((std::size_t)k * n, 0.0);
    for(int mode = 0; mode < k; mode++)
        for(int row = 0; row < n; row++)
            eigVector[(std::size_t)row + (std::size_t)mode * n] = vecs(row, mode);
    return true;
}
