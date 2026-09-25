#ifndef MICROGRID_H
#define MICROGRID_H

#include <stdint.h>
#include <stdio.h>
#include <float.h>
#include "MathPoint3D.h"
#include "MicroNode.h"
#include "MicroTetraedr.h"
#include "MicroFlags.h"
#include "FemMesh.h"

#define DEBUG_MODE
#define MIN_VOLUME_TETRAEDR (DBL_EPSILON)  /* Если объем тетраэдра меньше этой величины, то считать его вырожденным */
#define KOEF_SCALE_TETRAEDR (1.5)

class MicroGrid
{
public:
    MicroGrid();
    MicroGrid(const MathPoint3D& begin, const double& dx, const double& dy, const double& dz);
    void            setSizeGrid(const double& dx, const double& dy, const double& dz);
    double          getSizeX();
    double          getSizeY();
    double          getSizeZ();
    void            setBeginPoint(const MathPoint3D& point);
    MathPoint3D     getBeginPoint();
    MathPoint3D     getBeginModelPoint();
    MathPoint3D     getEndModelPoint();
    void            makeSuperStruct();
    bool            deleteSuperStruct();
    /* Load an externally generated conforming tetrahedral mesh (e.g. from Gmsh),
     * bypassing the built-in Delaunay generator. Builds nodes/tetrahedra, links
     * face-adjacency (surface convention nodes {i,(i+1)%4,(i+2)%4}), flags
     * boundary nodes NODE_IS_BOUNDARY, assigns serial numbers, and computes
     * boundary surfaces. After this the grid is ready for FEM assembly. */
    bool            loadTetraMesh(const FemMesh& mesh);
    bool            addNode(const MicroNode& node);
    bool            deleteOneTetraedr(ListMicroTetraedr::iterator &it);
    bool            linkNodes();
    void            clearLinkNodes();
    bool            setNumberNodesTetraedrs();
    bool            isAllDelone();
    bool            getIteratorNodes(std::list<MicroNode>::iterator& itBegin, std::list<MicroNode>::iterator& itEnd);
    bool            getIteratorTetraedrs(ListMicroTetraedr::iterator& itBegin, ListMicroTetraedr::iterator& itEnd);
    uint32_t        getNumNodes();
    uint32_t        getNumTetraedrs();
    double          volumeSuperStruct();
    void            clearTetraedrs();
    void            clear();
private:
    bool            tradeTetraedrs(ListMicroTetraedr::iterator tetraedr1, ListMicroTetraedr::iterator tetraedr2, ListMicroTetraedr::iterator *newIterator = 0);
    void            deleteTetraedr(ListMicroTetraedr::iterator &it);
    std::list<MicroNode> listNodes;
    ListMicroTetraedr listTetraedr;
    std::list<MicroTriangle> listTriangle;
    double          superStructDx;
    double          superStructDy;
    double          superStructDz;
    MathPoint3D     beginPoint;
    MathPoint3D     beginModelPoint, endModelPoint;

    bool            isClearTetraedrs;
    bool            isMakeSuperStruct;
    bool            isLinkNodes;
    bool            isNumNodesTetraedrs;
};


class MicroSurface
{
public:
    MicroSurface(ListMicroTetraedr::iterator itTetraedr, uint32_t numSurface):
        _tetraedr(itTetraedr), _numSurface(numSurface) {}
    ListMicroTetraedr::iterator&    tetraedr() {return _tetraedr;}
    uint32_t&                       numSurface() {return _numSurface;}
private:
    ListMicroTetraedr::iterator     _tetraedr;
    uint32_t                        _numSurface;
};

#endif // MICROGRID_H
