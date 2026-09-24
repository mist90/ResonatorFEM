#include "MicroNode.h"
#include "MicroTetraedr.h"
#include <stdlib.h>

/* пустой контейнер - для обнуления итераторов */
static std::list<MicroNode> emptyContainer;

MicroNode::MicroNode()
{
    _point = MathPoint3D(0, 0, 0);
    _serialNumber = 0;
    _neighbourNodes = 0;
    _areaFacesDirihle = 0;
    nNeighbourNodes = 0;
    nReservNeighbourNodes = 0;
    _flags = 0;
}

MicroNode::MicroNode(const MathPoint3D& point)
{
    _point = point;
    _serialNumber = 0;
    _neighbourNodes = 0;
    _areaFacesDirihle = 0;
    nNeighbourNodes = 0;
    nReservNeighbourNodes = 0;
    _flags = 0;
}

MicroNode::~MicroNode()
{
    clear();
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
    if((_flags & flags) == flags) return true;
    else return false;
}

uint32_t MicroNode::getSerialNumber()
{
    return _serialNumber;
}

bool MicroNode::addNeighbourNode(std::list<MicroNode>::iterator node)
{
    uint32_t i;
    for(i=0; i<nNeighbourNodes; i++)
        if(_neighbourNodes[i] == node) return false;
    if(nNeighbourNodes + 1 > nReservNeighbourNodes)
    {
        _neighbourNodes = (std::list<MicroNode>::iterator*)realloc(_neighbourNodes, (nReservNeighbourNodes + STEP_ALLOC) * sizeof(std::list<MicroNode>::iterator));
        _areaFacesDirihle = (double*)realloc(_areaFacesDirihle, (nReservNeighbourNodes + STEP_ALLOC) * sizeof(double));
        nReservNeighbourNodes += STEP_ALLOC;
    }
    _neighbourNodes[nNeighbourNodes] = node;
    _areaFacesDirihle[nNeighbourNodes] = 0.0;
    nNeighbourNodes++;
    return true;
}

void MicroNode::deleteNeighbourNode(uint32_t index)
{
    uint32_t i;
#ifdef DEBUG_MODE
    if(index >= nNeighbourNodes) exit(-1);
#endif
    for(i=index; i<nNeighbourNodes-1; i++)
    {
        _neighbourNodes[i] = _neighbourNodes[i+1];
        _areaFacesDirihle[i] = _areaFacesDirihle[i+1];
    }
    nNeighbourNodes--;
    nReservNeighbourNodes++;
}

bool MicroNode::isNeighbourNode(std::list<MicroNode>::iterator node)
{
    uint32_t i;
    for(i=0; i<numNeighbourNodes(); i++)
        if(neighbourNodes(i) == node) return true;
    return false;
}

std::list<MicroNode>::iterator& MicroNode::neighbourNodes(uint32_t index)
{
#ifdef DEBUG_MODE
    if(index >= nNeighbourNodes) exit(-1);
#endif
    return _neighbourNodes[index];
}

double MicroNode::areaFacesDirihle(uint32_t index)
{
#ifdef DEBUG_MODE
    if(index >= nNeighbourNodes) exit(-1);
#endif
    return _areaFacesDirihle[index];
}

double MicroNode::volumeDirihle()
{
    MathVector3D vector;
    uint32_t i;
    double volume = 0.0;
    for(i=0; i<nNeighbourNodes; i++)
    {
        vector = MathVector3D(point(), _neighbourNodes[i]->point());
        volume += (_areaFacesDirihle[i]*vector.lenght()/2.0);
    }
    volume = volume/3.0;
    return volume;
}

uint32_t MicroNode::numNeighbourNodes()
{
    return nNeighbourNodes;
}

void MicroNode::clearNeighbourNodes()
{
    if(_neighbourNodes)
    {
        free(_neighbourNodes);
        free(_areaFacesDirihle);
        _neighbourNodes = 0;
        _areaFacesDirihle = 0;
        nNeighbourNodes = 0;
    }
}

void MicroNode::clear()
{
    if(nNeighbourNodes)
    {
        free(_neighbourNodes);
        free(_areaFacesDirihle);
    }
    _neighbourNodes = 0;
    _areaFacesDirihle = 0;
    _flags = 0;
    nNeighbourNodes = 0;
}

std::list<MicroNode>::iterator EmptyNode()
{
    return emptyContainer.end();
}



//EOF
