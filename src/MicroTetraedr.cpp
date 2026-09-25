#include "MicroTetraedr.h"
#include "MathMatrix.h"
#include <stdlib.h>
#include <stdio.h>

/* пустой контейнер - для обнуления итераторов */
static ListMicroTetraedr emptyContainer;

/***********************************************************************************
 *                    class MicroEdge                                              *
 ***********************************************************************************/

MicroEdge::MicroEdge()
{
    _node1 = EmptyNode();
    _node2 = EmptyNode();
}

MicroEdge::MicroEdge(const std::list<MicroNode>::iterator &node1, const std::list<MicroNode>::iterator &node2)
{
    _node1 = node1;
    _node2 = node2;
}

void MicroEdge::setNode1(const std::list<MicroNode>::iterator &node1)
{
    _node1 = node1;
}

void MicroEdge::setNode2(const std::list<MicroNode>::iterator &node2)
{
    _node2 = node2;
}

std::list<MicroNode>::iterator MicroEdge::getNode1()
{
    return _node1;
}

std::list<MicroNode>::iterator MicroEdge::getNode2()
{
    return _node2;
}

bool MicroEdge::operator ==(const MicroEdge& edge)
{
    if(_node1 == edge._node1 && _node2 == edge._node2) return true;
    if(_node2 == edge._node1 && _node1 == edge._node2) return true;
    return false;
}

bool MicroEdge::isReverse(const MicroEdge &edge)
{
    if(_node1 == edge._node2 && _node2 == edge._node1) return true;
    else return false;
}

void MicroEdge::reverseEdge()
{
    std::list<MicroNode>::iterator itTmp;

    itTmp = _node1;
    _node1 = _node2;
    _node2 = itTmp;
}

double MicroEdge::lenEdge()
{
    return LenPoints(_node1->point(), _node2->point());
}

bool MicroEdge::isFlags(uint32_t flags)
{
    if(_node1->isFlags(flags) && _node2->isFlags(flags)) return true;
    else return false;
}


/***********************************************************************************
 *                    class MicroTetraedr                                            *
 ***********************************************************************************/

MicroTetraedr::MicroTetraedr()
{
    uint32_t i;
    for(i=0; i<4; i++)
    {
        _nodes[i] = EmptyNode();
        _tetraedrs[i] = EmptyIterator();
        _flags[i] = false;
    }
    _serialNumber = 0;
}

std::list<MicroNode>::iterator& MicroTetraedr::nodes(uint32_t i)
{
    return _nodes[i % 4];
}

ListMicroTetraedr::iterator& MicroTetraedr::tetraedrs(uint32_t i)
{
    return _tetraedrs[i % 4];
}

void MicroTetraedr::calculateBoundariesSurface()
{
    uint32_t i;
    for(i=0; i<4; i++)
        if(tetraedrs(i) == EmptyIterator())
        {
            if(nodes(i)->isFlags(NODE_IS_BOUNDARY) && nodes(i + 1)->isFlags(NODE_IS_BOUNDARY) && nodes(i + 2)->isFlags(NODE_IS_BOUNDARY))
                addFlags(i, SURFACE_IS_BOUNDARY);
        }
        else clearFlags(i, SURFACE_IS_BOUNDARY);
}

bool MicroTetraedr::isBoundarySurface(uint32_t index)
{
    return isFlags(index, SURFACE_IS_BOUNDARY);
}

bool MicroTetraedr::isFlags(uint32_t numSurface, uint32_t flags)
{
    numSurface = numSurface%4;
    if((_flags[numSurface] & flags) == flags) return true;
    else return false;
}

void MicroTetraedr::addFlags(uint32_t numSurface, uint32_t flags)
{
    numSurface = numSurface%4;
    _flags[numSurface] |= flags;
}

void MicroTetraedr::clearFlags(uint32_t numSurface, uint32_t flags)
{
    numSurface = numSurface%4;
    _flags[numSurface] &= (~flags);
}

double MicroTetraedr::volume()
{
    return VolumeTetraedr(nodes(0)->point(), nodes(1)->point(), nodes(2)->point(), nodes(3)->point());
}

void MicroTetraedr::setSerialNumber(uint32_t number)
{
    _serialNumber = number;
}

uint32_t MicroTetraedr::getSerialNumber()
{
    return _serialNumber;
}

/* Методы для реализации МКЭ */

uint32_t MicroTetraedr::getNumBeginEdge(uint32_t numEdge)
{
    numEdge = numEdge%6;
    switch(numEdge)
    {
    case 0:
        return 0;
        break;
    case 1:
        return 0;
        break;
    case 2:
        return 0;
        break;
    case 3:
        return 1;
        break;
    case 4:
        return 3;
        break;
    case 5:
        return 2;
        break;
    }
    return 0;
}

uint32_t MicroTetraedr::getNumEndEdge(uint32_t numEdge)
{
    numEdge = numEdge%6;
    switch(numEdge)
    {
    case 0:
        return 1;
        break;
    case 1:
        return 2;
        break;
    case 2:
        return 3;
        break;
    case 3:
        return 2;
        break;
    case 4:
        return 1;
        break;
    case 5:
        return 3;
        break;
    }
    return 0;
}

MicroEdge MicroTetraedr::getEdge(uint32_t numEdge)
{
    MicroEdge edge;

    numEdge = numEdge%6;
    edge.setNode1(nodes(getNumBeginEdge(numEdge)));
    edge.setNode2(nodes(getNumEndEdge(numEdge)));
    return edge;
}

MathVector3D MicroTetraedr::getValueBasisFunc(const MathPoint3D &point, uint32_t numEdge)
{
    MathMatrix<double> matrix(4, 4);
    MathVector3D vector1, vector2, ret;
    double a[2], b[2], c[2], d[2];  /* коэффициенты барицентрических функций */
    uint32_t i;

    for(i=0; i<4; i++)
    {
        matrix.element(i, 0) = nodes(i)->point().getX();
        matrix.element(i, 1) = nodes(i)->point().getY();
        matrix.element(i, 2) = nodes(i)->point().getZ();
        matrix.element(i, 3) = 1.0;
    }
    matrix = matrix.inverseMatrix();
    a[0] = matrix.element(3, getNumBeginEdge(numEdge));
    a[1] = matrix.element(3, getNumEndEdge(numEdge));
    b[0] = matrix.element(0, getNumBeginEdge(numEdge));
    b[1] = matrix.element(0, getNumEndEdge(numEdge));
    c[0] = matrix.element(1, getNumBeginEdge(numEdge));
    c[1] = matrix.element(1, getNumEndEdge(numEdge));
    d[0] = matrix.element(2, getNumBeginEdge(numEdge));
    d[1] = matrix.element(2, getNumEndEdge(numEdge));
    matrix.clear();
    vector1 = MathVector3D(b[1], c[1], d[1]);
    vector1 = vector1*(a[0] + b[0]*((MathPoint3D)point).getX() + c[0]*((MathPoint3D)point).getY() + d[0]*((MathPoint3D)point).getZ());
    vector2 = MathVector3D(b[0], c[0], d[0]);
    vector2 = vector2*(a[1] + b[1]*((MathPoint3D)point).getX() + c[1]*((MathPoint3D)point).getY() + d[1]*((MathPoint3D)point).getZ());
    ret = (vector1 - vector2)*getEdge(numEdge).lenEdge();
    ret.setBegin(point);
    return ret;
}

double VolumeTetraedr(MathPoint3D &point1, MathPoint3D &point2, MathPoint3D &point3, MathPoint3D &point4)
{
    MathMatrix<double> matrix;
    double volume;
    matrix.setSize(3, 3);
    matrix.element(0, 0) = point2.getX() - point1.getX();
    matrix.element(0, 1) = point3.getX() - point1.getX();
    matrix.element(0, 2) = point4.getX() - point1.getX();
    matrix.element(1, 0) = point2.getY() - point1.getY();
    matrix.element(1, 1) = point3.getY() - point1.getY();
    matrix.element(1, 2) = point4.getY() - point1.getY();
    matrix.element(2, 0) = point2.getZ() - point1.getZ();
    matrix.element(2, 1) = point3.getZ() - point1.getZ();
    matrix.element(2, 2) = point4.getZ() - point1.getZ();
    volume = fabs(matrix.determinantNoCopyMatrix()/6.0);
    return volume;
}
ListMicroTetraedr::ListMicroTetraedr(): std::list<MicroTetraedr>() { }

ListMicroTetraedr::iterator EmptyIterator()
{
    return emptyContainer.end();
}

//EOF
