#ifndef NODE_H
#define NODE_H
#include "MathPoint3D.h"
#include <stdint.h>
#include <list>

class MicroGrid;

/* A mesh vertex: coordinates, flags (boundary/metal), and a serial number. */
class MicroNode
{
public:
    MicroNode();
    MicroNode(const MathPoint3D& point);
    void            setPoint(const MathPoint3D& point);
    MathPoint3D&    point();
    void            addFlags(uint32_t flags);
    void            clearFlags(uint32_t flags);
    bool            isFlags(uint32_t flags);
    uint32_t        getSerialNumber();

private:
    MathPoint3D     _point;
    uint32_t        _flags;
    uint32_t        _serialNumber;
    friend class MicroGrid;   /* assigns _serialNumber in setNumberNodesTetraedrs */
};

std::list<MicroNode>::iterator EmptyNode();

#endif // NODE_H
