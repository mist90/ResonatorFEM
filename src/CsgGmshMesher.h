#ifndef CSGGMSHMESHER_H
#define CSGGMSHMESHER_H

#include <string>
#include "FemMesh.h"

/* Parametric resonator: a cavity (cylinder or box) with an optional metal core
 * (concentric cylinder, box, or elliptical-cross-section spiral) subtracted from
 * it. The core is centred on the cavity axis (z through the cavity centre). */
struct ResonatorSpec
{
    enum Cavity { CAVITY_CYLINDER = 0, CAVITY_BOX = 1 };
    enum Core   { CORE_NONE = 0, CORE_CYLINDER = 1, CORE_BOX = 2, CORE_SPIRAL = 3 };

    int    cavity   = CAVITY_CYLINDER;
    double cavRadius = 1.0, cavHeight = 2.0;      /* cylinder cavity            */
    double cavA = 1.0, cavB = 0.8, cavD = 1.5;    /* box cavity (x,y,z extents) */

    int    core = CORE_NONE;
    double coreRadius = 0.3, coreHeight = 2.0;    /* cylinder core              */
    double coreA = 0.3, coreB = 0.3, coreD = 1.0; /* box core (x,y,z extents)   */
    double helixR = 0.5, pitch = 0.5, turns = 2.0, ellA = 0.12, ellB = 0.12; /* spiral */

    double meshSize = 0.2;
};

/* Geometry + meshing backend built on Gmsh's OpenCASCADE (occ) kernel.
 *
 * Replaces the original hand-written Delaunay triangulator. Geometry is defined
 * with CSG primitives and boolean operations; Gmsh produces a conforming
 * tetrahedral mesh whose exterior faces become PEC walls in the FEM solver.
 *
 * Each call is self-contained: it (re)initializes Gmsh, builds the model, meshes
 * it, extracts a FemMesh, and clears the model. Not thread-safe (Gmsh holds
 * global state); drive it from a single thread. */
class CsgGmshMesher
{
public:
    /* Milestone 1: empty cylindrical cavity, axis along +z, base at origin.
     *   radius, height  - cavity dimensions (same length units used everywhere)
     *   meshSize        - target tetrahedron edge length
     * Returns true on success; on failure returns false and, if `err` != nullptr,
     * fills it with a message. */
    static bool cylindricalCavity(double radius, double height, double meshSize,
                                  FemMesh& out, std::string* err = nullptr);

    /* Milestone 1 (box variant): rectangular cavity spanning [0,a]x[0,b]x[0,d]. */
    static bool rectangularCavity(double a, double b, double d, double meshSize,
                                  FemMesh& out, std::string* err = nullptr);

    /* Milestone 2: coaxial cavity = outer cylinder MINUS a concentric inner
     * cylinder (the metal core), both along +z from the origin, full height.
     * The meshed region is the annulus between them; its inner and outer walls
     * and both end caps all become PEC. rInner < rOuter. */
    static bool coaxialCavity(double rInner, double rOuter, double height,
                              double meshSize, FemMesh& out, std::string* err = nullptr);

    /* Milestone 3: cylindrical cavity MINUS a helical (spiral) core with an
     * elliptical cross-section. The core is an ellipse (semi-axis ellA radial,
     * ellB axial) swept helically at radius helixR about the z-axis, `turns`
     * full turns with axial `pitch` per turn, centered in the cavity height.
     * The removed helix surface becomes PEC. Requires helixR+ellA < cavRadius
     * and pitch*turns + 2*ellB < cavHeight. */
    static bool spiralInCylindricalCavity(double cavRadius, double cavHeight,
                                          double helixR, double pitch, double turns,
                                          double ellA, double ellB, double meshSize,
                                          FemMesh& out, std::string* err = nullptr);

    /* Unified parametric builder: cavity {cyl, box} x core {none, cyl, box,
     * spiral}. Builds both solids, subtracts the core, meshes the remainder.
     * This is what the GUI drives. */
    static bool buildResonator(const ResonatorSpec& spec, FemMesh& out,
                               std::string* err = nullptr);
};

#endif // CSGGMSHMESHER_H
