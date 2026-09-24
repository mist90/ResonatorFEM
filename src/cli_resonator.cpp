/* Headless validation harness (milestone 1).
 *
 * Meshes an empty cavity with the Gmsh/OCC backend, runs the reused edge-element
 * FEM resonator eigen-solve, and prints the resonant frequencies alongside the
 * analytic reference for a cylindrical cavity.
 *
 * Usage:
 *   resonator_cli [cyl|box] [dim...] [meshSize]
 *     cyl radius height meshSize     (default: 1 2 0.6)
 *     box a b d meshSize
 */
#include <QCoreApplication>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "CsgGmshMesher.h"
#include "MicroEngine.h"

static const double C_LIGHT = 299792458.0;         /* m/s */

/* n-th positive zero of the Bessel function J_0 (for TM_0m0 modes). */
static const double J0_ZERO_1 = 2.4048255576957727;

static double k2_to_MHz(double k2)
{
    if (k2 < 0) return 0.0;
    return C_LIGHT * std::sqrt(k2) / (2.0 * M_PI) / 1.0e6;
}

int main(int argc, char** argv)
{
    QCoreApplication app(argc, argv);

    std::string shape = (argc > 1) ? argv[1] : "cyl";
    FemMesh mesh;
    std::string err;
    double analyticTM010_MHz = 0.0;

    if (shape == "spiral") {
        double cavR = 1.0, cavH = 2.0;
        double helixR = (argc > 2) ? std::atof(argv[2]) : 0.5;
        double pitch  = (argc > 3) ? std::atof(argv[3]) : 0.5;
        double turns  = (argc > 4) ? std::atof(argv[4]) : 2.0;
        double ellA   = (argc > 5) ? std::atof(argv[5]) : 0.12;
        double ellB   = (argc > 6) ? std::atof(argv[6]) : 0.12;
        double m      = (argc > 7) ? std::atof(argv[7]) : 0.18;
        std::printf("Spiral core in cyl cavity: cavR=%.2f cavH=%.2f helixR=%.2f pitch=%.2f "
                    "turns=%.1f ellipse=%.2fx%.2f  meshSize=%.3f\n",
                    cavR, cavH, helixR, pitch, turns, ellA, ellB, m);
        if (!CsgGmshMesher::spiralInCylindricalCavity(cavR, cavH, helixR, pitch, turns,
                                                      ellA, ellB, m, mesh, &err)) {
            std::printf("meshing FAILED: %s\n", err.c_str());
            return 1;
        }
        /* Rough helical quarter-wave estimate (wire length quarter-wavelength);
         * real helical resonators deviate due to inter-turn loading. */
        double wireLen = turns * std::sqrt(std::pow(2.0 * M_PI * helixR, 2.0) + pitch * pitch);
        analyticTM010_MHz = C_LIGHT / (4.0 * wireLen) / 1.0e6;
    } else if (shape == "coax") {
        double ri = (argc > 2) ? std::atof(argv[2]) : 0.3;
        double ro = (argc > 3) ? std::atof(argv[3]) : 1.0;
        double h  = (argc > 4) ? std::atof(argv[4]) : 2.0;
        double m  = (argc > 5) ? std::atof(argv[5]) : 0.35;
        std::printf("Coaxial cavity rInner=%.3f rOuter=%.3f height=%.3f  meshSize=%.3f\n",
                    ri, ro, h, m);
        if (!CsgGmshMesher::coaxialCavity(ri, ro, h, m, mesh, &err)) {
            std::printf("meshing FAILED: %s\n", err.c_str());
            return 1;
        }
        /* Shorted-TEM fundamental: f = c/(2h), independent of the radii. */
        analyticTM010_MHz = C_LIGHT / (2.0 * h) / 1.0e6;
    } else if (shape == "box") {
        double a = (argc > 2) ? std::atof(argv[2]) : 1.0;
        double b = (argc > 3) ? std::atof(argv[3]) : 0.8;
        double d = (argc > 4) ? std::atof(argv[4]) : 1.5;
        double h = (argc > 5) ? std::atof(argv[5]) : 0.6;
        std::printf("Rectangular cavity a=%.3f b=%.3f d=%.3f  meshSize=%.3f\n", a, b, d, h);
        if (!CsgGmshMesher::rectangularCavity(a, b, d, h, mesh, &err)) {
            std::printf("meshing FAILED: %s\n", err.c_str());
            return 1;
        }
        /* TE101 (dominant for a box, if a>=b, d>=b): f=c/2*sqrt(1/a^2+1/d^2) */
        analyticTM010_MHz = C_LIGHT / 2.0 *
            std::sqrt(1.0 / (a * a) + 1.0 / (d * d)) / 1.0e6;
    } else {
        double r = (argc > 2) ? std::atof(argv[2]) : 1.0;
        double h = (argc > 3) ? std::atof(argv[3]) : 2.0;
        double m = (argc > 4) ? std::atof(argv[4]) : 0.6;
        std::printf("Cylindrical cavity radius=%.3f height=%.3f  meshSize=%.3f\n", r, h, m);
        if (!CsgGmshMesher::cylindricalCavity(r, h, m, mesh, &err)) {
            std::printf("meshing FAILED: %s\n", err.c_str());
            return 1;
        }
        /* TM010: f = c * chi01 / (2*pi*r) */
        analyticTM010_MHz = C_LIGHT * J0_ZERO_1 / (2.0 * M_PI * r) / 1.0e6;
    }
    std::printf("mesh: %zu nodes, %zu tets\n", mesh.nodes.size(), mesh.tets.size());

    MicroEngine engine;                 /* headless (no grapher) */
    engine.setResonatorMode(true);
    /* Sparse shift-invert targeting the expected physical range (k^2 of the
     * analytic/estimate frequency): isolates physical modes, skips the null
     * space, and scales to fine meshes. */
    double targetK2 = std::pow(2.0 * M_PI * analyticTM010_MHz * 1.0e6 / C_LIGHT, 2.0);
    double sigmaK2 = 0.9 * targetK2;   /* sit in the gap just below the target mode */
    engine.setResonatorSolveTarget(sigmaK2, 12);
    std::printf("shift-invert: target %.1f MHz (k^2=%.4f), shift sigma=%.4f\n",
                analyticTM010_MHz, targetK2, sigmaK2);
    if (!engine.generateFromMesh(mesh)) {
        std::printf("generateFromMesh FAILED\n");
        return 1;
    }
    std::printf("solving generalized eigenproblem (sparse shift-invert, Spectra)...\n");
    std::fflush(stdout);
    if (!engine.calculateSync()) {
        std::printf("eigen-solve FAILED\n");
        return 1;
    }

    std::vector<double> k2;
    if (!engine.getEighValues(k2) || k2.empty()) {
        std::printf("no eigenvalues returned\n");
        return 1;
    }
    std::sort(k2.begin(), k2.end());

    /* Count the near-zero gradient null-space, then print the full non-null
     * spectrum for direct comparison with the analytic mode table. */
    const double NULL_TOL = 1.0e-2;
    size_t nNull = 0;
    while (nNull < k2.size() && k2[nNull] < NULL_TOL) ++nNull;

    std::printf("\ntotal DOFs (free edges): %zu, null-space modes (k^2<%.0e): %zu\n",
                k2.size(), NULL_TOL, nNull);
    std::printf("\n  # |        k^2 |    f (MHz) | note\n");
    std::printf("----+------------+------------+-----------------------\n");
    double bestErr = 1e9;
    for (size_t i = nNull, shown = 0; i < k2.size() && shown < 20; ++i, ++shown) {
        double f = k2_to_MHz(k2[i]);
        double err = 100.0 * (f - analyticTM010_MHz) / analyticTM010_MHz;
        const char* note = "";
        if (std::fabs(err) < std::fabs(bestErr)) { bestErr = err; }
        if (std::fabs(err) < 5.0) note = "<- near analytic dominant";
        std::printf("%3zu | %10.5f | %10.3f | %s\n", shown, k2[i], f, note);
    }
    std::printf("\nanalytic dominant mode: %.3f MHz   (closest computed mode error %+.2f%%)\n",
                analyticTM010_MHz, bestErr);
    return 0;
}
