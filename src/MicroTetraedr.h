#ifndef MICROTETRAEDR_H
#define MICROTETRAEDR_H

#include "MathVector3D.h"
#include "MathVector3D.h"
#include "MicroNode.h"
#include "MicroFlags.h"
#include <stdint.h>
#include <list>
#include <vector>

/* An edge of a tetrahedron (an ordered pair of node iterators). */
class MicroEdge
{
public:
    MicroEdge();
    MicroEdge(const std::list<MicroNode>::iterator& node1, const std::list<MicroNode>::iterator& node2);
    void                            setNode1(const std::list<MicroNode>::iterator& node1);
    void                            setNode2(const std::list<MicroNode>::iterator& node2);
    std::list<MicroNode>::iterator  getNode1();
    std::list<MicroNode>::iterator  getNode2();
    /* True if the two edges share the same node pair (either orientation). */
    bool                            operator==(const MicroEdge& edge);
    /* True if the two edges are the same pair but opposite orientation. */
    bool                            isReverse(const MicroEdge& edge);
    double                          lenEdge();
private:
    std::list<MicroNode>::iterator _node1;
    std::list<MicroNode>::iterator _node2;
};

class MicroTetraedr;
/* Container of tetrahedra. */
class ListMicroTetraedr:public std::list<MicroTetraedr>
{
public:
    ListMicroTetraedr();
};

/* A tetrahedron: 4 node iterators, 4 neighbour links, per-surface flags, and the
 * edge-element (Whitney) FEM helpers. */
class MicroTetraedr
{
public:
    MicroTetraedr();
    std::list<MicroNode>::iterator& nodes(uint32_t i);
    ListMicroTetraedr::iterator&    tetraedrs(uint32_t i);
    /* Mark boundary surfaces (no neighbour); call after the tet is fully built. */
    void                            calculateBoundariesSurface();
    bool                            isBoundarySurface(uint32_t index);
    /* Per-surface flags. */
    bool                            isFlags(uint32_t numSurface, uint32_t flags);
    void                            addFlags(uint32_t numSurface, uint32_t flags);
    void                            clearFlags(uint32_t numSurface, uint32_t flags);
    double                          volume();
    void                            setSerialNumber(uint32_t number);
    uint32_t                        getSerialNumber();
    /* Edge enumeration (6 edges) for FEM assembly. */
    uint32_t                        getNumBeginEdge(uint32_t numEdge);
    uint32_t                        getNumEndEdge(uint32_t numEdge);
    MicroEdge                       getEdge(uint32_t numEdge);
    /* Value of Whitney basis function numEdge at a point. */
    MathVector3D                    getValueBasisFunc(const MathVector3D& point, uint32_t numEdge);

private:
    std::list<MicroNode>::iterator      _nodes[4];
    ListMicroTetraedr::iterator         _tetraedrs[4];
    uint32_t                            _flags[4];
    uint32_t                            _serialNumber;
};

/* An "empty" tetrahedron iterator sentinel (used for missing neighbours). */
ListMicroTetraedr::iterator EmptyIterator();

/* Signed-magnitude volume of a tetrahedron. */
double VolumeTetraedr(MathVector3D &point1, MathVector3D &point2, MathVector3D &point3, MathVector3D &point4);

#endif // MICROTETRAEDR_H
