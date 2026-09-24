#ifndef MATHEIGHVALUES
#define MATHEIGHVALUES
#include "MathMatrixSparse.h"
#include <vector>
#include <utility>
#include <cstdint>

#define KOEF_MAX_DIV    (100.0)
//#define DEBUG_MODE

/* A*X=l*B*X, dense LAPACK dsygv: full spectrum, small/well-conditioned problems. */
bool MathEighValVector(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB, std::vector<double> &eighValue, std::vector<double>& eigVector);

/* A*X=l*B*X, sparse shift-invert (Spectra): the `nev` eigenvalues nearest `sigma`.
 * Targets the physical mode range directly (skipping the near-zero edge-element
 * null space) and only factorizes (A - sigma*B), so it tolerates a singular mass
 * matrix B. Eigenvectors packed column-major: eigVector[row + mode*n]. */
bool MathEighValVectorShiftInvert(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB,
                                  double sigma, int nev,
                                  std::vector<double>& eighValue, std::vector<double>& eigVector);

/* As above, plus a mass-metric grad-div penalty A += penaltyS * (T G) D^-1 (Gᵀ T),
 * where G is the discrete gradient given by dofNodes[dof] = (nodeSerial1,nodeSerial2)
 * (UINT32_MAX entries are skipped) and D = diag(Gᵀ T G). The penalty is exactly
 * zero on physical modes (Gᵀ T e = 0) and lifts the gradient null-space modes to
 * ~penaltyS, so the physical modes become the smallest eigenvalues. */
bool MathEighValVectorShiftInvertGauged(MathMatrixSparse<double>& matrixA, MathMatrixSparse<double>& matrixB,
                                        const std::vector<std::pair<uint32_t, uint32_t> >& dofNodes,
                                        const std::vector<double>& dofLen,
                                        uint32_t nNodes, double penaltyS, double sigma, int nev,
                                        std::vector<double>& eighValue, std::vector<double>& eigVector);

#endif // MATHEIGHVALUES

