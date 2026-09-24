#include "CsgGmshMesher.h"

#include <gmsh.h>
#include <cstddef>
#include <exception>
#include <cmath>
#include <unordered_map>

namespace {

/* Configure common meshing options (quiet, target element size). */
void applyMeshOptions(double meshSize)
{
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::option::setNumber("Mesh.MeshSizeMin", meshSize * 0.5);
    gmsh::option::setNumber("Mesh.MeshSizeMax", meshSize);
    /* Delaunay 3D is robust and fast for these simple CSG solids. */
    gmsh::option::setNumber("Mesh.Algorithm3D", 1);
    /* Optimize element quality: sliver tetrahedra near curved (faceted) walls
     * perturb the edge-element gradient null space toward small nonzero values,
     * so we remove them. */
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.OptimizeThreshold", 0.5);
}

/* Extract the current (already-generated) 3D mesh into `out`.
 * Node tags from Gmsh are arbitrary/1-based; we remap them to dense 0-based
 * indices. Each volume entity contributes its tets with a distinct region tag
 * (0,1,2,... in entity order) so the FEM layer can assign per-region material. */
bool extractCurrentMesh(FemMesh& out, std::string* err)
{
    out.clear();

    std::vector<std::size_t> nodeTags;
    std::vector<double>      coords;
    std::vector<double>      paramCoords;
    gmsh::model::mesh::getNodes(nodeTags, coords, paramCoords);

    if (nodeTags.empty()) {
        if (err) *err = "mesh has no nodes";
        return false;
    }

    std::unordered_map<std::size_t, uint32_t> tagToIndex;
    tagToIndex.reserve(nodeTags.size() * 2);
    out.nodes.resize(nodeTags.size());
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        tagToIndex[nodeTags[i]] = static_cast<uint32_t>(i);
        out.nodes[i] = { coords[3 * i + 0], coords[3 * i + 1], coords[3 * i + 2] };
    }

    /* Type 4 == 4-node linear tetrahedron. Gather per volume entity so region
     * tags follow the CSG solids. */
    std::vector<std::pair<int, int> > volumes;
    gmsh::model::getEntities(volumes, 3);
    int region = 0;
    for (const auto& v : volumes) {
        std::vector<std::size_t> elemTags, elemNodeTags;
        gmsh::model::mesh::getElementsByType(4, elemTags, elemNodeTags, v.second);
        const std::size_t nTet = elemTags.size();
        for (std::size_t e = 0; e < nTet; ++e) {
            std::array<uint32_t, 4> tet;
            for (int k = 0; k < 4; ++k)
                tet[k] = tagToIndex.at(elemNodeTags[4 * e + k]);
            out.tets.push_back(tet);
            out.tetRegion.push_back(region);
        }
        ++region;
    }

    if (out.tets.empty()) {
        if (err) *err = "mesh has no tetrahedra";
        return false;
    }
    return true;
}

} // namespace

bool CsgGmshMesher::cylindricalCavity(double radius, double height, double meshSize,
                                      FemMesh& out, std::string* err)
{
    try {
        gmsh::initialize();
        gmsh::model::add("cyl_cavity");
        /* Cylinder: base center (0,0,0), axis (0,0,height), radius. */
        gmsh::model::occ::addCylinder(0, 0, 0, 0, 0, height, radius);
        gmsh::model::occ::synchronize();
        applyMeshOptions(meshSize);
        gmsh::model::mesh::generate(3);
        bool ok = extractCurrentMesh(out, err);
        gmsh::clear();
        gmsh::finalize();
        return ok;
    } catch (const std::exception& e) {
        if (err) *err = std::string("gmsh: ") + e.what();
        try { gmsh::finalize(); } catch (...) {}
        return false;
    }
}

bool CsgGmshMesher::rectangularCavity(double a, double b, double d, double meshSize,
                                      FemMesh& out, std::string* err)
{
    try {
        gmsh::initialize();
        gmsh::model::add("box_cavity");
        gmsh::model::occ::addBox(0, 0, 0, a, b, d);
        gmsh::model::occ::synchronize();
        applyMeshOptions(meshSize);
        gmsh::model::mesh::generate(3);
        bool ok = extractCurrentMesh(out, err);
        gmsh::clear();
        gmsh::finalize();
        return ok;
    } catch (const std::exception& e) {
        if (err) *err = std::string("gmsh: ") + e.what();
        try { gmsh::finalize(); } catch (...) {}
        return false;
    }
}

bool CsgGmshMesher::coaxialCavity(double rInner, double rOuter, double height,
                                  double meshSize, FemMesh& out, std::string* err)
{
    if (rInner <= 0.0 || rInner >= rOuter) {
        if (err) *err = "coaxialCavity: require 0 < rInner < rOuter";
        return false;
    }
    try {
        gmsh::initialize();
        gmsh::model::add("coax_cavity");
        int outer = gmsh::model::occ::addCylinder(0, 0, 0, 0, 0, height, rOuter);
        int inner = gmsh::model::occ::addCylinder(0, 0, 0, 0, 0, height, rInner);
        /* annulus = outer - inner (the metal core is removed from the domain). */
        gmsh::vectorpair outDimTags;
        std::vector<gmsh::vectorpair> outDimTagsMap;
        gmsh::model::occ::cut({ {3, outer} }, { {3, inner} }, outDimTags, outDimTagsMap);
        gmsh::model::occ::synchronize();
        applyMeshOptions(meshSize);
        gmsh::model::mesh::generate(3);
        bool ok = extractCurrentMesh(out, err);
        gmsh::clear();
        gmsh::finalize();
        return ok;
    } catch (const std::exception& e) {
        if (err) *err = std::string("gmsh: ") + e.what();
        try { gmsh::finalize(); } catch (...) {}
        return false;
    }
}

bool CsgGmshMesher::spiralInCylindricalCavity(double cavRadius, double cavHeight,
                                              double helixR, double pitch, double turns,
                                              double ellA, double ellB, double meshSize,
                                              FemMesh& out, std::string* err)
{
    const double height = pitch * turns;
    if (helixR + ellA >= cavRadius || helixR - ellA <= 0.0) {
        if (err) *err = "spiral: core does not fit radially (need 0 < helixR-ellA and helixR+ellA < cavRadius)";
        return false;
    }
    if (height + 2.0 * ellB >= cavHeight) {
        if (err) *err = "spiral: core too tall (need pitch*turns + 2*ellB < cavHeight)";
        return false;
    }
    try {
        gmsh::initialize();
        gmsh::model::add("spiral_cavity");

        int cavity = gmsh::model::occ::addCylinder(0, 0, 0, 0, 0, cavHeight, cavRadius);

        const double z0 = (cavHeight - height) / 2.0;
        const double totalAngle = turns * 2.0 * M_PI;

        /* Helix spine as an OCC B-spline through sampled points. */
        const int perTurn = 24;
        const int nPts = (int)(perTurn * turns) + 1;
        std::vector<int> pts;
        pts.reserve(nPts);
        for (int i = 0; i < nPts; ++i) {
            double t = totalAngle * (double)i / (double)(nPts - 1);
            double x = helixR * std::cos(t);
            double y = helixR * std::sin(t);
            double z = z0 + pitch * t / (2.0 * M_PI);
            pts.push_back(gmsh::model::occ::addPoint(x, y, z, meshSize));
        }
        int spline = gmsh::model::occ::addBSpline(pts);
        int wire = gmsh::model::occ::addWire({ spline });

        /* Elliptical profile at the helix start, normal = helix tangent there.
         * Tangent at t=0 is (0, helixR, pitch/2pi); radial x-axis (1,0,0) lies in
         * the profile plane. Major semi-axis ellA is radial, minor ellB axial. */
        double tz = pitch / (2.0 * M_PI);
        double tn = std::sqrt(helixR * helixR + tz * tz);
        std::vector<double> zAxis = { 0.0, helixR / tn, tz / tn };  /* = unit tangent */
        std::vector<double> xAxis = { 1.0, 0.0, 0.0 };
        int ell = gmsh::model::occ::addEllipse(helixR, 0.0, z0, ellA, ellB, -1,
                                               0.0, 2.0 * M_PI, zAxis, xAxis);
        int loop = gmsh::model::occ::addCurveLoop({ ell });
        int face = gmsh::model::occ::addPlaneSurface({ loop });

        /* Sweep the elliptical face along the helix -> helical elliptical solid. */
        gmsh::vectorpair swept;
        gmsh::model::occ::addPipe({ {2, face} }, wire, swept, "CorrectedFrenet");
        int core = -1;
        for (const auto& dt : swept)
            if (dt.first == 3) { core = dt.second; break; }
        if (core < 0) {
            if (err) *err = "spiral: pipe sweep did not produce a solid";
            gmsh::clear(); gmsh::finalize();
            return false;
        }

        gmsh::vectorpair outDimTags;
        std::vector<gmsh::vectorpair> outDimTagsMap;
        gmsh::model::occ::cut({ {3, cavity} }, { {3, core} }, outDimTags, outDimTagsMap);
        gmsh::model::occ::synchronize();
        applyMeshOptions(meshSize);
        gmsh::model::mesh::generate(3);
        bool ok = extractCurrentMesh(out, err);
        gmsh::clear();
        gmsh::finalize();
        return ok;
    } catch (const std::exception& e) {
        if (err) *err = std::string("gmsh: ") + e.what();
        try { gmsh::finalize(); } catch (...) {}
        return false;
    }
}
