#ifndef CSGGMSHMESHER_H
#define CSGGMSHMESHER_H

#include <string>
#include "FemMesh.h"

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
};

#endif // CSGGMSHMESHER_H
