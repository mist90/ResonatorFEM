#ifndef MATHEIGHVALUES
#define MATHEIGHVALUES
#include "MathMatrixSparse.h"
#include <vector>

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

#endif // MATHEIGHVALUES

