#include "MathEighValues.h"
#include <vector>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <random>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
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

bool MathEighValVector(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                       std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.rows() != matrixA.cols()) return false;
    if(matrixB.rows() != matrixB.cols()) return false;
    if(matrixA.rows() != matrixB.rows()) return false;

    const int n = (int)matrixA.rows();
    if(n <= 0) return false;

    /* Dense column-major copies (Eigen is column-major); dsygv overwrites A with
     * the eigenvectors and B with the Cholesky factor. */
    Eigen::MatrixXd A = Eigen::MatrixXd(matrixA);
    Eigen::MatrixXd B = Eigen::MatrixXd(matrixB);

    eighValue.assign(n, 0.0);
    const int itype = 1;          /* A x = lambda B x            */
    const char jobz = 'V';        /* eigenvalues and eigenvectors */
    const char uplo = 'U';
    int lda = n, ldb = n, info = 0, lwork = -1;

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

    /* A now holds the eigenvectors as columns (column-major): eigVector[row+mode*n]. */
    eigVector.assign(A.data(), A.data() + (std::size_t)n * n);
    return true;
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

/* Truncated Rayleigh-Ritz: the `k` smallest Ritz pairs of the small dense pencil
 * (gramA, gramB) on an m-dimensional subspace, dropping directions where gramB is
 * near-singular (the [X W P] basis becomes linearly dependent as LOBPCG
 * converges). Returns Ritz values `theta` and coefficients `C` (m x k') with
 * k' = min(k, numeric rank); Ritz vectors are V*C and are B-orthonormal. */
static bool truncatedRR(const Eigen::MatrixXd& gramA, const Eigen::MatrixXd& gramB,
                        int k, Eigen::VectorXd& theta, Eigen::MatrixXd& C)
{
    const int m = (int)gramA.rows();
    if(m < 1) return false;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esB(gramB);
    if(esB.info() != Eigen::Success) return false;
    const Eigen::VectorXd& bev = esB.eigenvalues();     /* ascending */
    double bmax = bev(m - 1);
    if(bmax <= 0.0) return false;
    double thr = bmax * 1e-12;
    std::vector<int> keep;
    for(int i = 0; i < m; i++) if(bev(i) > thr) keep.push_back(i);
    int r = (int)keep.size();
    if(r < 1) return false;
    /* T maps the kept gramB-eigenbasis to a gramB-orthonormal basis. */
    Eigen::MatrixXd T(m, r);
    for(int j = 0; j < r; j++)
        T.col(j) = esB.eigenvectors().col(keep[j]) / std::sqrt(bev(keep[j]));
    Eigen::MatrixXd Ared = T.transpose() * gramA * T;
    Ared = 0.5 * (Ared + Ared.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esA(Ared);
    if(esA.info() != Eigen::Success) return false;
    int kk = std::min(k, r);
    theta = esA.eigenvalues().head(kk);
    C = T * esA.eigenvectors().leftCols(kk);            /* m x kk */
    return true;
}

/* Block LOBPCG: the `nev` smallest eigenpairs of A x = lambda B x with A, B SPD.
 * Iterative and matrix-free apart from an incomplete-Cholesky preconditioner on A
 * (falls back to Jacobi) — no direct factorization, so memory stays ~O(nnz(A))
 * instead of the shift-invert factor fill-in. Meant for the grad-div-penalized
 * operator Aeff, which is SPD and far more elliptic than raw curl-curl. */
static bool lobpcgSmallest(const Eigen::SparseMatrix<double>& A,
                           const Eigen::SparseMatrix<double>& B,
                           int nev, int maxIter, double tol,
                           Eigen::VectorXd& vals, Eigen::MatrixXd& vecs)
{
    const int n = (int)A.rows();
    int k = nev;
    if(k > n / 3) k = n / 3;
    if(k < 1) return false;

    /* Preconditioner K ~ A^-1: incomplete Cholesky, else Jacobi. */
    Eigen::IncompleteCholesky<double> ic;
    ic.compute(A);
    bool haveIC = (ic.info() == Eigen::Success);
    Eigen::VectorXd invdiag;
    if(!haveIC) {
        invdiag.resize(n);
        for(int i = 0; i < n; i++) { double d = A.coeff(i, i); invdiag[i] = (d > 1e-30) ? 1.0 / d : 1.0; }
        std::fprintf(stderr, "lobpcg: IncompleteCholesky failed -> Jacobi preconditioner\n");
    }
    auto precond = [&](const Eigen::MatrixXd& Rm) -> Eigen::MatrixXd {
        Eigen::MatrixXd Wm(n, Rm.cols());
        if(haveIC) for(int c = 0; c < Rm.cols(); c++) Wm.col(c) = ic.solve(Rm.col(c));
        else       for(int c = 0; c < Rm.cols(); c++) Wm.col(c) = invdiag.cwiseProduct(Rm.col(c));
        return Wm;
    };

    /* Deterministic random start (mt19937 fixed seed: reproducible). */
    std::mt19937 gen(1234567u);
    std::normal_distribution<double> nd(0.0, 1.0);
    Eigen::MatrixXd X(n, k);
    for(int j = 0; j < k; j++) for(int i = 0; i < n; i++) X(i, j) = nd(gen);

    Eigen::MatrixXd AX = A * X, BX = B * X;
    Eigen::VectorXd theta; Eigen::MatrixXd C;
    {
        Eigen::MatrixXd gA = X.transpose() * AX, gB = X.transpose() * BX;
        gA = 0.5 * (gA + gA.transpose()); gB = 0.5 * (gB + gB.transpose());
        if(!truncatedRR(gA, gB, k, theta, C)) return false;
        X = X * C; AX = AX * C; BX = BX * C;
    }
    int kc = (int)theta.size();
    Eigen::VectorXd Lam = theta;
    Eigen::MatrixXd R = AX - BX * Lam.asDiagonal();
    Eigen::MatrixXd P;   /* search direction block (0 cols initially) */

    int iter = 0;
    double maxrel = 1.0;
    for(; iter < maxIter; ++iter) {
        maxrel = 0.0;
        for(int c = 0; c < kc; c++)
            maxrel = std::max(maxrel, R.col(c).norm() / (AX.col(c).norm() + 1e-30));
        if(maxrel < tol) break;

        Eigen::MatrixXd W = precond(R);
        int p = (int)P.cols();
        int m = 2 * kc + p;
        Eigen::MatrixXd V(n, m);
        V.leftCols(kc) = X;
        V.middleCols(kc, kc) = W;
        if(p > 0) V.rightCols(p) = P;
        Eigen::MatrixXd AV = A * V, BV = B * V;
        Eigen::MatrixXd gA = V.transpose() * AV, gB = V.transpose() * BV;
        gA = 0.5 * (gA + gA.transpose()); gB = 0.5 * (gB + gB.transpose());
        if(!truncatedRR(gA, gB, kc, theta, C)) break;   /* keep last good X */
        int knew = (int)theta.size();
        Eigen::MatrixXd Cwp = C.bottomRows(m - kc);      /* [W P] contribution */
        P  = V.rightCols(m - kc) * Cwp;                  /* new search direction */
        X  = V * C; AX = AV * C; BX = BV * C;
        Lam = theta; kc = knew;
        R = AX - BX * Lam.asDiagonal();
    }
    std::fprintf(stderr, "lobpcg: %d iters, %d modes, max rel-resid %.2e (%s), precond=%s\n",
                 iter, kc, maxrel, maxrel < tol ? "converged" : "NOT converged",
                 haveIC ? "IC" : "Jacobi");
    vals = Lam;
    vecs = X;
    return kc >= 1;
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

bool MathEighValVectorShiftInvert(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                                  double sigma, int nev,
                                  std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.rows() != matrixA.cols()) return false;
    if(matrixB.rows() != matrixB.cols()) return false;
    if(matrixA.rows() != matrixB.rows()) return false;
    const int n = (int)matrixA.rows();
    if(n <= 2) return false;

    Eigen::VectorXd vals;
    Eigen::MatrixXd vecs;
    if(!shiftInvertCore(matrixA, matrixB, sigma, nev, vals, vecs)) return false;
    packResults(vals, vecs, n, eighValue, eigVector);
    return true;
}

bool MathEighValVectorShiftInvertGauged(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                                        const std::vector<std::pair<uint32_t, uint32_t> >& dofNodes,
                                        const std::vector<double>& dofLen,
                                        const std::vector<char>& interiorNode,
                                        uint32_t nNodes, double penaltyS, double sigma, int nev,
                                        std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.rows() != matrixA.cols()) return false;
    if(matrixB.rows() != matrixB.cols()) return false;
    if(matrixA.rows() != matrixB.rows()) return false;
    const int n = (int)matrixA.rows();
    if(n <= 2 || (int)dofNodes.size() != n || (int)dofLen.size() != n || nNodes == 0) return false;
    if((uint32_t)interiorNode.size() != nNodes) return false;

    const Eigen::SparseMatrix<double>& Ae = matrixA;   /* stiffness S */
    const Eigen::SparseMatrix<double>& Be = matrixB;   /* mass T      */

    /* Discrete gradient G (n edges x nNodes). For this len*Whitney basis the DOF
     * coefficient of grad(phi) on edge (tail a -> head b) is (phi_b - phi_a)/len,
     * so entries are -/+ 1/len. Crucially, only INTERIOR nodes are columns: the
     * valid gauge freedom is phi varying on interior nodes and constant on each
     * PEC conductor. The gradient of a boundary-node scalar leaks onto the
     * excluded PEC edges and is NOT a null vector of the reduced stiffness, so
     * including boundary-node columns spoils G (this was the key fix). With the
     * interior-only restriction G is (essentially) the null space of S, so the
     * penalty preserves the physical modes and lifts only the gradient modes. */
    std::vector<Eigen::Triplet<double> > gtrip;
    gtrip.reserve((std::size_t)n * 2);
    for(int d = 0; d < n; d++) {
        uint32_t a = dofNodes[d].first, b = dofNodes[d].second;
        if(a == UINT32_MAX || dofLen[d] <= 0.0) continue;
        double w = 1.0 / dofLen[d];
        if(interiorNode[a]) gtrip.emplace_back(d, (int)a, -w);
        if(interiorNode[b]) gtrip.emplace_back(d, (int)b, +w);
    }
    Eigen::SparseMatrix<double> G((int)n, (int)nNodes);
    G.setFromTriplets(gtrip.begin(), gtrip.end());
    /* Gauge-quality report: ||S*G||/||S|| should be ~0 (G in null(S)). It is
     * ~0.01-0.07 for cavities; larger on coarse meshes of complex curved
     * boundaries (e.g. the spiral) and shrinks under refinement. */
    { Eigen::SparseMatrix<double> SG = Ae * G;
      std::fprintf(stderr, "gauge: ||S*G||/||S|| = %.3e\n", SG.norm() / (Ae.norm() + 1e-30)); }

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
    /* Solver: default is Spectra shift-invert (direct factorization, robust but
     * memory-heavy on fine 3D meshes). RESONATOR_SOLVER=lobpcg selects the
     * iterative low-memory LOBPCG; if it fails it falls back to shift-invert. */
    const char* solver = std::getenv("RESONATOR_SOLVER");
    if(solver && std::strcmp(solver, "lobpcg") == 0) {
        int    maxIter = 800;
        double tol     = 1e-4;   /* eigenvalue error ~ tol^2, so ~1e-8 — plenty */
        if(const char* mi = std::getenv("LOBPCG_ITERS")) maxIter = std::atoi(mi);
        if(const char* tl = std::getenv("LOBPCG_TOL"))   tol     = std::atof(tl);
        if(lobpcgSmallest(Aeff, Be, nev, maxIter, tol, vals, vecs)) {
            packResults(vals, vecs, n, eighValue, eigVector);
            return true;
        }
        std::fprintf(stderr, "lobpcg failed -> falling back to shift-invert\n");
    }
    if(!shiftInvertCore(Aeff, Be, sigma, nev, vals, vecs)) return false;
    packResults(vals, vecs, n, eighValue, eigVector);
    return true;
}
