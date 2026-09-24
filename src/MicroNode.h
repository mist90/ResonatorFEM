#ifndef NODE_H
#define NODE_H
#include "MathPoint3D.h"
#include <stdint.h>
#include <list>

#define STEP_ALLOC (6)
#define DEBUG_MODE

class MicroTetraedr;
class MicroGrid;

class MicroNode
{
public:
    MicroNode();
    MicroNode(const MathPoint3D& point);
    ~MicroNode();
    void            setPoint(const MathPoint3D& point);
    MathPoint3D&    point();
    void            addFlags(uint32_t flags);
    void            clearFlags(uint32_t flags);
    bool            isFlags(uint32_t flags);
    uint32_t        getSerialNumber();
    bool            addNeighbourNode(std::list<MicroNode>::iterator node);
    void            deleteNeighbourNode(uint32_t index);
    bool            isNeighbourNode(std::list<MicroNode>::iterator node);
    std::list<MicroNode>::iterator&     neighbourNodes(uint32_t index);
    double          areaFacesDirihle(uint32_t index);
    double          volumeDirihle();
    uint32_t        numNeighbourNodes();
    void            clearNeighbourNodes();
    void            clear();

private:
    MathPoint3D     _point;
    uint32_t        _flags;
    uint32_t        _serialNumber;
    std::list<MicroNode>::iterator*     _neighbourNodes;    /* не используется std::list для экономии памяти */
    double*         _areaFacesDirihle;
    uint32_t        nNeighbourNodes;
    uint32_t        nReservNeighbourNodes;
    friend class MicroGrid;
};

std::list<MicroNode>::iterator EmptyNode();


#endif // NODE_H
