#include "MicroNode.h"

/* Empty container, for a null node iterator sentinel. */
static std::list<MicroNode> emptyContainer;

MicroNode::MicroNode()
{
    _point = MathPoint3D(0, 0, 0);
    _flags = 0;
    _serialNumber = 0;
}

MicroNode::MicroNode(const MathPoint3D& point)
{
    _point = point;
    _flags = 0;
    _serialNumber = 0;
}

void MicroNode::setPoint(const MathPoint3D& point)
{
    _point = point;
}

MathPoint3D& MicroNode::point()
{
    return _point;
}

void MicroNode::addFlags(uint32_t flags)
{
    _flags |= flags;
}

void MicroNode::clearFlags(uint32_t flags)
{
    _flags &= (~flags);
}

bool MicroNode::isFlags(uint32_t flags)
{
    return (_flags & flags) == flags;
}

uint32_t MicroNode::getSerialNumber()
{
    return _serialNumber;
}

std::list<MicroNode>::iterator EmptyNode()
{
    return emptyContainer.end();
}
