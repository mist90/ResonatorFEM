#ifndef MICROGRID_H
#define MICROGRID_H

#include <stdint.h>
#include "MicroNode.h"
#include "MicroTetraedr.h"
#include "MicroFlags.h"
#include "FemMesh.h"

/* Container of the FEM mesh (nodes + tetrahedra). The geometry is produced
 * externally by Gmsh and loaded via loadTetraMesh; this class no longer
 * generates a mesh itself. */
class MicroGrid
{
public:
    MicroGrid();

    /* Load an externally generated conforming tetrahedral mesh (from Gmsh):
     * build nodes/tetrahedra, link face-adjacency (surface convention nodes
     * {i,(i+1)%4,(i+2)%4}), flag boundary nodes NODE_IS_BOUNDARY, assign serial
     * numbers, and compute boundary surfaces. Ready for FEM assembly after. */
    bool            loadTetraMesh(const FemMesh& mesh);

    bool            setNumberNodesTetraedrs();
    bool            getIteratorNodes(std::list<MicroNode>::iterator& itBegin, std::list<MicroNode>::iterator& itEnd);
    bool            getIteratorTetraedrs(ListMicroTetraedr::iterator& itBegin, ListMicroTetraedr::iterator& itEnd);
    uint32_t        getNumNodes();
    uint32_t        getNumTetraedrs();
    void            clear();

private:
    std::list<MicroNode>  listNodes;
    ListMicroTetraedr     listTetraedr;
    bool                  isClearTetraedrs;
    bool                  isMakeSuperStruct;
    bool                  isLinkNodes;
    bool                  isNumNodesTetraedrs;
};

#endif // MICROGRID_H
