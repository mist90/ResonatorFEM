#include "MathEighValues.h"
#include <mpi.h>
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_parcsr_ls.h>   /* HYPRE_AME* are declared only in the internal header */
#include <vector>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>

/* hypre AME (Auxiliary-space Maxwell Eigensolver): block LOBPCG preconditioned
 * with AMS, the auxiliary-space AMG for H(curl). AME keeps the iterates discretely
 * divergence-free (M-orthogonal to range(G)), so it works on the raw, singular
 * curl-curl operator and returns only the physical modes — no grad-div penalty.
 *
 * hypre expects the standard lowest-order Nedelec DOFs (edge circulations) and
 * the +-1 discrete gradient G. The engine's basis is length-scaled Whitney,
 * N_e = len_e * W_e, so its coefficients are x_e = y_e / len_e with y the
 * standard DOFs. With D = diag(len) the pencil becomes
 *     (D^-1 A D^-1) y = lambda (D^-1 B D^-1) y,   x = D^-1 y,
 * which has the same eigenvalues.
 *
 * Boundary conditions follow hypre's (and MFEM's) convention: the system covers
 * every mesh edge; PEC edges are kept as identity rows of A and zero rows of M.
 * AME detects those rows and removes the boundary vertices from G itself. */

namespace {

typedef std::vector<std::pair<uint32_t, uint32_t> > EdgeList;

/* Build a (rows x cols) ParCSR matrix from CSR arrays through the IJ interface. */
HYPRE_IJMatrix makeIJMatrix(HYPRE_BigInt rows, HYPRE_BigInt cols,
                            const std::vector<HYPRE_Int>& rowNnz,
                            std::vector<HYPRE_BigInt>& colIdx,
                            std::vector<HYPRE_Real>& vals)
{
    HYPRE_IJMatrix ij;
    HYPRE_IJMatrixCreate(MPI_COMM_SELF, 0, rows - 1, 0, cols - 1, &ij);
    HYPRE_IJMatrixSetObjectType(ij, HYPRE_PARCSR);
    std::vector<HYPRE_Int> sizes(rowNnz);
    HYPRE_IJMatrixSetRowSizes(ij, sizes.data());
    HYPRE_IJMatrixInitialize(ij);
    std::vector<HYPRE_Int>    ncols(rowNnz);
    std::vector<HYPRE_BigInt> rowIds((std::size_t)rows);
    for(HYPRE_BigInt i = 0; i < rows; i++) rowIds[(std::size_t)i] = i;
    HYPRE_IJMatrixSetValues(ij, (HYPRE_Int)rows, ncols.data(), rowIds.data(),
                            colIdx.data(), vals.data());
    HYPRE_IJMatrixAssemble(ij);
    return ij;
}

HYPRE_IJVector makeIJVector(const std::vector<HYPRE_Real>& v)
{
    HYPRE_IJVector ij;
    const HYPRE_BigInt n = (HYPRE_BigInt)v.size();
    HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, n - 1, &ij);
    HYPRE_IJVectorSetObjectType(ij, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(ij);
    std::vector<HYPRE_BigInt> idx((std::size_t)n);
    for(HYPRE_BigInt i = 0; i < n; i++) idx[(std::size_t)i] = i;
    std::vector<HYPRE_Real> vals(v);
    HYPRE_IJVectorSetValues(ij, (HYPRE_Int)n, idx.data(), vals.data());
    HYPRE_IJVectorAssemble(ij);
    return ij;
}

/* Full-edge operator: the scaled free-DOF block S (rows 0..n-1), then `nPec`
 * eliminated rows holding `pecDiag` on the diagonal. */
HYPRE_IJMatrix makeEdgeOperator(const Eigen::SparseMatrix<double, Eigen::RowMajor>& S,
                                int nPec, double pecDiag)
{
    const int n = (int)S.rows();
    const HYPRE_BigInt N = (HYPRE_BigInt)n + nPec;
    std::vector<HYPRE_Int>    rowNnz((std::size_t)N);
    std::vector<HYPRE_BigInt> colIdx;
    std::vector<HYPRE_Real>   vals;
    colIdx.reserve((std::size_t)S.nonZeros() + nPec);
    vals.reserve((std::size_t)S.nonZeros() + nPec);
    for(int r = 0; r < n; r++) {
        HYPRE_Int cnt = 0;
        for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(S, r); it; ++it) {
            colIdx.push_back(it.col());
            vals.push_back(it.value());
            cnt++;
        }
        rowNnz[(std::size_t)r] = cnt;
    }
    for(int k = 0; k < nPec; k++) {
        colIdx.push_back((HYPRE_BigInt)n + k);
        vals.push_back(pecDiag);
        rowNnz[(std::size_t)n + k] = 1;
    }
    return makeIJMatrix(N, N, rowNnz, colIdx, vals);
}

int envInt(const char* name, int def)
{
    const char* v = std::getenv(name);
    return v ? std::atoi(v) : def;
}

double envDouble(const char* name, double def)
{
    const char* v = std::getenv(name);
    return v ? std::atof(v) : def;
}

} // namespace

bool MathEighValVectorAME(const Eigen::SparseMatrix<double>& matrixA, const Eigen::SparseMatrix<double>& matrixB,
                          const EdgeList& dofNodes,
                          const std::vector<double>& dofLen,
                          const EdgeList& pecEdges,
                          const std::vector<std::array<double, 3> >& nodeCoords,
                          int nev,
                          std::vector<double>& eighValue, std::vector<double>& eigVector)
{
    if(matrixA.rows() != matrixA.cols() || matrixA.rows() != matrixB.rows() ||
       matrixB.rows() != matrixB.cols()) return false;
    const int n = (int)matrixA.rows();
    const int nPec = (int)pecEdges.size();
    const int nv = (int)nodeCoords.size();
    if(n <= 2 || nev < 1 || nev >= n || nv == 0) return false;
    if((int)dofNodes.size() != n || (int)dofLen.size() != n) return false;
    for(int d = 0; d < n; d++)
        if(dofLen[d] <= 0.0 || dofNodes[d].first >= (uint32_t)nv || dofNodes[d].second >= (uint32_t)nv)
            return false;

    int mpiInit = 0;
    MPI_Initialized(&mpiInit);
    if(!mpiInit) {
        std::fprintf(stderr, "ame: MPI/hypre not initialized (no EigenSolverRuntime)\n");
        return false;
    }

    /* Rescale to standard Nedelec DOFs (see the file comment). */
    Eigen::VectorXd dinv(n);
    for(int d = 0; d < n; d++) dinv[d] = 1.0 / dofLen[d];
    Eigen::SparseMatrix<double, Eigen::RowMajor> As = dinv.asDiagonal() * matrixA * dinv.asDiagonal();
    Eigen::SparseMatrix<double, Eigen::RowMajor> Bs = dinv.asDiagonal() * matrixB * dinv.asDiagonal();
    As.makeCompressed();
    Bs.makeCompressed();

    HYPRE_ClearAllErrors();
    HYPRE_IJMatrix ijA = makeEdgeOperator(As, nPec, 1.0);   /* PEC rows: identity */
    HYPRE_IJMatrix ijM = makeEdgeOperator(Bs, nPec, 0.0);   /* PEC rows: zero     */

    /* Discrete gradient over all edges and vertices: row e = -1 at its first
     * node, +1 at its second (the DOF's reference orientation). */
    const HYPRE_BigInt N = (HYPRE_BigInt)n + nPec;
    HYPRE_IJMatrix ijG;
    {
        std::vector<HYPRE_Int>    rowNnz((std::size_t)N, 2);
        std::vector<HYPRE_BigInt> colIdx;
        std::vector<HYPRE_Real>   vals;
        colIdx.reserve((std::size_t)N * 2);
        vals.reserve((std::size_t)N * 2);
        for(HYPRE_BigInt e = 0; e < N; e++) {
            const std::pair<uint32_t, uint32_t>& en =
                (e < n) ? dofNodes[(std::size_t)e] : pecEdges[(std::size_t)(e - n)];
            colIdx.push_back(en.first);  vals.push_back(-1.0);
            colIdx.push_back(en.second); vals.push_back(+1.0);
        }
        ijG = makeIJMatrix(N, nv, rowNnz, colIdx, vals);
    }

    std::vector<HYPRE_Real> cx(nv), cy(nv), cz(nv), zeros((std::size_t)N, 0.0);
    for(int i = 0; i < nv; i++) { cx[i] = nodeCoords[i][0]; cy[i] = nodeCoords[i][1]; cz[i] = nodeCoords[i][2]; }
    HYPRE_IJVector ijX = makeIJVector(cx), ijY = makeIJVector(cy), ijZ = makeIJVector(cz);
    HYPRE_IJVector ijB = makeIJVector(zeros), ijX0 = makeIJVector(zeros);

    HYPRE_ParCSRMatrix A, M, G;
    HYPRE_ParVector x, y, z, b, x0;
    HYPRE_IJMatrixGetObject(ijA, (void**)&A);
    HYPRE_IJMatrixGetObject(ijM, (void**)&M);
    HYPRE_IJMatrixGetObject(ijG, (void**)&G);
    HYPRE_IJVectorGetObject(ijX, (void**)&x);
    HYPRE_IJVectorGetObject(ijY, (void**)&y);
    HYPRE_IJVectorGetObject(ijZ, (void**)&z);
    HYPRE_IJVectorGetObject(ijB, (void**)&b);
    HYPRE_IJVectorGetObject(ijX0, (void**)&x0);

    /* AMS as a single-V-cycle preconditioner. Beta = 0 (no mass term in A:
     * singular curl-curl). Parameters follow hypre's ams_driver / MFEM defaults. */
    HYPRE_Solver ams;
    HYPRE_AMSCreate(&ams);
    HYPRE_AMSSetDimension(ams, 3);
    HYPRE_AMSSetMaxIter(ams, 1);
    HYPRE_AMSSetTol(ams, 0.0);
    HYPRE_AMSSetCycleType(ams, 13);
    HYPRE_AMSSetPrintLevel(ams, 0);
    HYPRE_AMSSetDiscreteGradient(ams, G);
    HYPRE_AMSSetCoordinateVectors(ams, x, y, z);
    HYPRE_AMSSetBetaPoissonMatrix(ams, NULL);
    HYPRE_AMSSetSmoothingOptions(ams, 2, 1, 1.0, 1.0);
    HYPRE_AMSSetAlphaAMGOptions(ams, 10, 1, 8, 0.25, 6, 4);
    HYPRE_AMSSetBetaAMGOptions(ams, 10, 1, 8, 0.25, 6, 4);
    HYPRE_AMSSetAlphaAMGCoarseRelaxType(ams, 8);
    HYPRE_AMSSetBetaAMGCoarseRelaxType(ams, 8);
    HYPRE_AMSSetup(ams, A, b, x0);

    const int    maxIter = envInt("AME_ITERS", 200);
    const double tol     = envDouble("AME_TOL", 1e-6);
    HYPRE_Solver ame;
    HYPRE_AMECreate(&ame);
    HYPRE_AMESetAMSSolver(ame, ams);
    HYPRE_AMESetMassMatrix(ame, M);
    HYPRE_AMESetBlockSize(ame, nev);
    HYPRE_AMESetMaxIter(ame, maxIter);
    HYPRE_AMESetTol(ame, tol);
    HYPRE_AMESetPrintLevel(ame, 0);
    HYPRE_AMESetup(ame);
    HYPRE_AMESolve(ame);
    /* The inner PCG projection may flag HYPRE_ERROR_CONV without harm; the
     * result is validated below by its residual instead. */
    HYPRE_ClearAllErrors();

    HYPRE_Real*      evals = NULL;
    HYPRE_ParVector* evecs = NULL;
    HYPRE_AMEGetEigenvalues(ame, &evals);
    HYPRE_AMEGetEigenvectors(ame, &evecs);

    /* Back to the engine's basis: x = D^-1 y (free DOFs only). */
    Eigen::VectorXd vals(nev);
    Eigen::MatrixXd vecs(n, nev);
    for(int k = 0; k < nev; k++) {
        vals[k] = evals[k];
        const HYPRE_Real* data = hypre_VectorData(hypre_ParVectorLocalVector((hypre_ParVector*)evecs[k]));
        for(int d = 0; d < n; d++) vecs(d, k) = data[d] * dinv[d];
        HYPRE_ParVectorDestroy(evecs[k]);
    }
    hypre_TFree(evecs, HYPRE_MEMORY_HOST);
    hypre_TFree(evals, HYPRE_MEMORY_HOST);

    HYPRE_AMEDestroy(ame);
    HYPRE_AMSDestroy(ams);
    HYPRE_IJVectorDestroy(ijX0); HYPRE_IJVectorDestroy(ijB);
    HYPRE_IJVectorDestroy(ijZ);  HYPRE_IJVectorDestroy(ijY); HYPRE_IJVectorDestroy(ijX);
    HYPRE_IJMatrixDestroy(ijG);  HYPRE_IJMatrixDestroy(ijM); HYPRE_IJMatrixDestroy(ijA);

    /* Validate with the original pencil: rel-resid = ||A x - l B x|| / (s ||B x||),
     * s = max(|l|, largest |l|). Scaling by the block's eigenvalue scale keeps a
     * genuine k^2 ~ 0 mode (electrostatic mode of a floating conductor) from
     * reading as unconverged. */
    double lScale = 0.0;
    for(int k = 0; k < nev; k++) lScale = std::max(lScale, std::fabs(vals[k]));
    double maxRes = 0.0;
    bool finite = true;
    for(int k = 0; k < nev; k++) {
        Eigen::VectorXd v  = vecs.col(k);
        Eigen::VectorXd Bv = matrixB * v;
        double denom = std::max(std::fabs(vals[k]), lScale) * Bv.norm();
        double res = (denom > 0.0) ? (matrixA * v - vals[k] * Bv).norm() / denom : INFINITY;
        if(!std::isfinite(vals[k]) || !std::isfinite(res)) finite = false;
        maxRes = std::max(maxRes, res);
    }
    std::fprintf(stderr, "ame: %d modes, %d PEC edges, max rel-resid %.2e%s\n",
                 nev, nPec, maxRes, (maxRes > 10.0 * tol) ? " (NOT converged)" : "");
    if(!finite || maxRes > 1e-2) return false;

    /* Ascending order, eigVector[row + mode*n] layout. */
    std::vector<int> order(nev);
    for(int k = 0; k < nev; k++) order[k] = k;
    std::sort(order.begin(), order.end(), [&vals](int a, int b) { return vals[a] < vals[b]; });
    eighValue.resize(nev);
    eigVector.assign((std::size_t)nev * n, 0.0);
    for(int k = 0; k < nev; k++) {
        eighValue[k] = vals[order[k]];
        for(int d = 0; d < n; d++)
            eigVector[(std::size_t)d + (std::size_t)k * n] = vecs(d, order[k]);
    }
    return true;
}

static bool g_ownsMpi = false;

EigenSolverRuntime::EigenSolverRuntime(int* argc, char*** argv)
{
    int inited = 0;
    MPI_Initialized(&inited);
    if(!inited) {
        /* SERIALIZED: the GUI solves on a worker thread, one solve at a time. */
        int provided = 0;
        MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);
        g_ownsMpi = true;
    }
    HYPRE_Init();
}

EigenSolverRuntime::~EigenSolverRuntime()
{
    HYPRE_Finalize();
    int finalized = 0;
    MPI_Finalized(&finalized);
    if(g_ownsMpi && !finalized) MPI_Finalize();
}
