#ifndef FEMMESH_H
#define FEMMESH_H

#include <array>
#include <cstdint>
#include <vector>

/* Plain, Qt-free tetrahedral mesh interchange format.
 * Produced by the geometry/mesh backend (CsgGmshMesher) and consumed by the
 * FEM mesh adapter (MicroGrid::loadTetraMesh). Node indices in `tets` are
 * 0-based offsets into `nodes`. `tetRegion` (optional, same length as `tets`)
 * carries a material/region tag per tetrahedron for per-element permittivity. */
struct FemMesh
{
    std::vector<std::array<double, 3> >   nodes;      /* xyz coordinates          */
    std::vector<std::array<uint32_t, 4> > tets;       /* 4 node indices per tet   */
    std::vector<int>                      tetRegion;  /* region tag per tet (opt) */

    void clear()
    {
        nodes.clear();
        tets.clear();
        tetRegion.clear();
    }
    bool empty() const { return tets.empty(); }
};

#endif // FEMMESH_H
