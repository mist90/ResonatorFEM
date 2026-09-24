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
    _pointCenterSphere = MathPoint3D(0.0, 0.0, 0.0);
    _radSphereQuad = 0.0;
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

void MicroTetraedr::clearTetraedrs()
{
    uint32_t i;
    for(i=0; i<4; i++)
    {
        _nodes[i] = EmptyNode();
        _tetraedrs[i] = EmptyIterator();
        _flags[i] = false;
    }
}

bool MicroTetraedr::isInTetraedr(const MathPoint3D &point)
{
    uint32_t i;

    for(i=0; i<4; i++)
        if(!isInTetraedr(point, i))
            return false;
    return true;
}

bool MicroTetraedr::isInTetraedr(const MathPoint3D &point, uint32_t numSurface)
{
    uint32_t i;
    double m, n, p;
    double dt1, dt2;

    MathPoint3D *points[4];
    for(i=0; i<4; i++) if(_nodes[i] == EmptyNode()) return false;
    for(i=0; i<4; i++) points[i] = &_nodes[(numSurface + i)%4]->point();

    m = (points[1]->getY() - points[0]->getY())*(points[2]->getZ() - points[0]->getZ()) -
            (points[1]->getZ() - points[0]->getZ())*(points[2]->getY() - points[0]->getY());
    n = -(points[1]->getX() - points[0]->getX())*(points[2]->getZ() - points[0]->getZ()) +
            (points[1]->getZ() - points[0]->getZ())*(points[2]->getX() - points[0]->getX());
    p = (points[1]->getX() - points[0]->getX())*(points[2]->getY() - points[0]->getY()) -
            (points[1]->getY() - points[0]->getY())*(points[2]->getX() - points[0]->getX());
    if(m*m + n*n + p*p == 0) return false;
    dt1 = m*(points[3]->getX() - points[0]->getX()) + n*(points[3]->getY() - points[0]->getY()) + p*(points[3]->getZ() - points[0]->getZ());
    dt2 = m*(((MathPoint3D)point).getX() - points[0]->getX()) + n*(((MathPoint3D)point).getY() - points[0]->getY()) + p*(((MathPoint3D)point).getZ() - points[0]->getZ());
    if(dt1*dt2 < 0) return false;

    return true;
}

bool MicroTetraedr::isPointInSurface(const MathPoint3D &point, uint32_t numSurface)
{
    uint32_t i;
    for(i=0; i<3; i++)
    {
        if(nodes(numSurface + i) != EmptyNode())
            if(nodes(numSurface + i)->point() == point) return true;
    }
    return false;
}

bool MicroTetraedr::isEdgeInSurface(const MathPoint3D &point1, const MathPoint3D &point2, uint32_t numSurface)
{
    if(!isPointInSurface(point1, numSurface)) return false;
    if(!isPointInSurface(point2, numSurface)) return false;
    return true;
}

bool MicroTetraedr::isNodeTetraedr(const MathPoint3D &point)
{
    uint32_t i;
    for(i=0; i<4; i++)
        if(_nodes[i] != EmptyNode())
        {
            if(_nodes[i]->point() == point) return true;
        }

    return false;
}

bool MicroTetraedr::isNode3Tetraedr(const MathPoint3D &point1, const MathPoint3D &point2, const MathPoint3D &point3)
{
    if(!isNodeTetraedr(point1)) return false;
    if(!isNodeTetraedr(point2)) return false;
    if(!isNodeTetraedr(point3)) return false;
    return true;
}

bool MicroTetraedr::getNumSurface(const MathPoint3D &point1, const MathPoint3D &point2, const MathPoint3D &point3, uint32_t &getNum)
{
    std::list<uint32_t> listSurfaces;
    std::list<uint32_t>::iterator listSurfacesIterator;
    uint32_t i;

    for(i=0; i<4; i++) listSurfaces.push_back(i);
    for(listSurfacesIterator=listSurfaces.begin(); listSurfacesIterator!=listSurfaces.end(); listSurfacesIterator++)
        if(!isPointInSurface(point1, *listSurfacesIterator))
        {
            listSurfaces.erase(listSurfacesIterator);
            break;
        }
    for(listSurfacesIterator=listSurfaces.begin(); listSurfacesIterator!=listSurfaces.end(); listSurfacesIterator++)
        if(!isPointInSurface(point2, *listSurfacesIterator))
        {
            listSurfaces.erase(listSurfacesIterator);
            break;
        }
    for(listSurfacesIterator=listSurfaces.begin(); listSurfacesIterator!=listSurfaces.end(); listSurfacesIterator++)
        if(!isPointInSurface(point3, *listSurfacesIterator))
        {
            listSurfaces.erase(listSurfacesIterator);
            break;
        }
    if(listSurfaces.size() != 1) return false;
    else getNum = *listSurfaces.begin();

    return true;
}

bool MicroTetraedr::getNumSurface(ListMicroTetraedr::iterator it, uint32_t &getNum)
{
    uint32_t i;

    for(i=0; i<4; i++)
        if(tetraedrs(i) == it)
        {
            getNum = i;
            return true;
        }
    return false;
}

uint32_t MicroTetraedr::getCountEmptyNeighbours()
{
    uint32_t i, count = 0;
    for(i=0; i<4; i++) if(tetraedrs(i) == EmptyIterator()) count++;
    return count;
}

void MicroTetraedr::makeCircumSphere()
{
#ifndef CALC_SPHERE_METHOD_SLAU
    uint32_t i;
    double a, c, Dx, Dy, Dz;
    double sumQuad[4];
    MathMatrix<double> matrix = MathMatrix<double>(4, 4);

    for(i=0; i<4; i++) if(_nodes[i] == EmptyNode()) return;
    for(i=0; i<4; i++) sumQuad[i] = nodes(i)->point().getX()*nodes(i)->point().getX()+
                                    nodes(i)->point().getY()*nodes(i)->point().getY()+
                                    nodes(i)->point().getZ()*nodes(i)->point().getZ();
    for(i=0; i<4; i++)
    {
        matrix.element(0, i) = nodes(i)->point().getX();
        matrix.element(1, i) = nodes(i)->point().getY();
        matrix.element(2, i) = nodes(i)->point().getZ();
        matrix.element(3, i) = 1;
    }
    a = matrix.determinantNoCopyMatrix();

    for(i=0; i<4; i++)
    {
        matrix.element(0, i) = sumQuad[i];
        matrix.element(1, i) = nodes(i)->point().getX();
        matrix.element(2, i) = nodes(i)->point().getY();
        matrix.element(3, i) = nodes(i)->point().getZ();
    }
    c = matrix.determinantNoCopyMatrix();

    for(i=0; i<4; i++)
    {
        matrix.element(0, i) = sumQuad[i];
        matrix.element(1, i) = nodes(i)->point().getY();
        matrix.element(2, i) = nodes(i)->point().getZ();
        matrix.element(3, i) = 1;
    }
    Dx = matrix.determinantNoCopyMatrix();
    for(i=0; i<4; i++)
    {
        matrix.element(0, i) = sumQuad[i];
        matrix.element(1, i) = nodes(i)->point().getX();
        matrix.element(2, i) = nodes(i)->point().getZ();
        matrix.element(3, i) = 1;
    }
    Dy = -matrix.determinantNoCopyMatrix();
    for(i=0; i<4; i++)
    {
        matrix.element(0, i) = sumQuad[i];
        matrix.element(1, i) = nodes(i)->point().getX();
        matrix.element(2, i) = nodes(i)->point().getY();
        matrix.element(3, i) = 1;
    }
    Dz = matrix.determinantNoCopyMatrix();
    _pointCenterSphere.setX(Dx/2.0/a);
    _pointCenterSphere.setY(Dy/2.0/a);
    _pointCenterSphere.setZ(Dz/2.0/a);
    _radSphereQuad = (Dx*Dx + Dy*Dy + Dz*Dz - 4.0*a*c)/4.0/(a*a);

#else

    MathMatrix<double> matrix = MathMatrix<double>(4, 3);
    std::vector<double> roots;
    double radQuad[4];
    uint32_t i;

    for(i=0; i<3; i++)
    {
        matrix.element(0, i) = nodes(i + 1)->point().getX() - nodes(i)->point().getX();
        matrix.element(1, i) = nodes(i + 1)->point().getY() - nodes(i)->point().getY();
        matrix.element(2, i) = nodes(i + 1)->point().getZ() - nodes(i)->point().getZ();
        matrix.element(3, i) = 0.5*(nodes(i + 1)->point().getX()*nodes(i + 1)->point().getX() - nodes(i)->point().getX()*nodes(i)->point().getX() +
                                    nodes(i + 1)->point().getY()*nodes(i + 1)->point().getY() - nodes(i)->point().getY()*nodes(i)->point().getY() +
                                    nodes(i + 1)->point().getZ()*nodes(i + 1)->point().getZ() - nodes(i)->point().getZ()*nodes(i)->point().getZ());
    }

    matrix.solveNoCopyMatrix(roots);
    _pointCenterSphere.setX(roots[0]);
    _pointCenterSphere.setY(roots[1]);
    _pointCenterSphere.setZ(roots[2]);

    for(i=0; i<4; i++)
    {
        radQuad[i] = (nodes(i)->point().getX() - _pointCenterSphere.getX())*(nodes(i)->point().getX() - _pointCenterSphere.getX()) +
                (nodes(i)->point().getY() - _pointCenterSphere.getY())*(nodes(i)->point().getY() - _pointCenterSphere.getY()) +
                (nodes(i)->point().getZ() - _pointCenterSphere.getZ())*(nodes(i)->point().getZ() - _pointCenterSphere.getZ());
    }
    _radSphereQuad = (radQuad[0] + radQuad[1] + radQuad[2] + radQuad[3])/4.0;
#endif

}

MathPoint3D MicroTetraedr::pointCenterSphere()
{
    return _pointCenterSphere;
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

bool MicroTetraedr::isInSphere(const MathPoint3D &point)
{
    double lenPointsQuad;

    lenPointsQuad = (((MathPoint3D)point).getX() - _pointCenterSphere.getX())*(((MathPoint3D)point).getX() - _pointCenterSphere.getX()) +
            (((MathPoint3D)point).getY() - _pointCenterSphere.getY())*(((MathPoint3D)point).getY() - _pointCenterSphere.getY()) +
            (((MathPoint3D)point).getZ() - _pointCenterSphere.getZ())*(((MathPoint3D)point).getZ() - _pointCenterSphere.getZ());
    if(lenPointsQuad < _radSphereQuad*(1.0 - EPSILON_RAD_QUAD))
        return true;
    else
        return false;
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

uint32_t MicroTetraedr::getNumEdgeFromSurface(uint32_t numSurface, uint32_t numEdge)
{
    numSurface = numSurface%4;
    numEdge = numEdge%3;
    switch (numSurface) {
    case 0:
        switch (numEdge) {
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return 3;
        }

    case 1:
        switch (numEdge) {
        case 0:
            return 0;
        case 1:
            return 2;
        case 2:
            return 4;
        }

    case 2:
        switch (numEdge) {
        case 0:
            return 1;
        case 1:
            return 2;
        case 2:
            return 5;
        }

    case 3:
        switch (numEdge) {
        case 0:
            return 3;
        case 1:
            return 4;
        case 2:
            return 5;
        }

    }
    return 0;
}

uint32_t MicroTetraedr::getNumPointFromSurface(uint32_t numSurface, uint32_t numPoint)
{
    numSurface = numSurface%4;
    numPoint = numPoint%3;
    switch (numSurface) {
    case 0:
        switch (numPoint) {
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return 2;
        }

    case 1:
        switch (numPoint) {
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return 3;
        }

    case 2:
        switch (numPoint) {
        case 0:
            return 0;
        case 1:
            return 2;
        case 2:
            return 3;
        }

    case 3:
        switch (numPoint) {
        case 0:
            return 1;
        case 1:
            return 2;
        case 2:
            return 3;
        }

    }
    return 0;
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

MathVector3D MicroTetraedr::getNormVector(uint32_t numSurface)
{
    MathVector3D ret;
    MathVector3D vector1 = MathVector3D(nodes(numSurface)->point(), nodes(numSurface + 1)->point());
    MathVector3D vector2 = MathVector3D(nodes(numSurface)->point(), nodes(numSurface + 2)->point());

    ret = vector1^vector2;
    ret = ret.unitVector();
    ret.setBegin((nodes(numSurface)->point() + nodes(numSurface + 1)->point() + nodes(numSurface + 2)->point())/3.0);
    if(isInTetraedr(ret.getEnd(), numSurface)) ret = -ret;
    return ret;
}


/************************************************************************************
 *                           class ListMicroTetraedr                                *
 ************************************************************************************/
ListMicroTetraedr::ListMicroTetraedr(): std::list<MicroTetraedr>() { }

bool ListMicroTetraedr::findTetraedr(MathPoint3D &point, ListMicroTetraedr::iterator &it)
{
    ListMicroTetraedr::iterator itTmp;
    for(itTmp=this->begin(); itTmp!=this->end(); itTmp++)
    {
        if(itTmp->isInTetraedr(point))
        {
            it = itTmp;
            return true;
        }
    }
    return false;    /* Ошибка: не удалось найти тетраэдр */
}

ListMicroTetraedr::iterator EmptyIterator()
{
    return emptyContainer.end();
}

/***********************************************************************************
 *                    class MicroTriangle                                          *
 ***********************************************************************************/
MicroTriangle::MicroTriangle(std::list<MicroNode>::iterator node1, std::list<MicroNode>::iterator node2, std::list<MicroNode>::iterator node3)
{
    _nodes[0] = node1;
    _nodes[1] = node2;
    _nodes[2] = node3;
}

std::list<MicroNode>::iterator MicroTriangle::node(uint32_t index)
{
    return _nodes[index % 3];
}

bool MicroTriangle::isNodeTriangle(MathPoint3D &point)
{
    uint32_t i;
    for(i=0; i<3; i++)
        if(_nodes[i] != EmptyNode())
            if(_nodes[i]->point() == point) return true;
    return false;
}

bool MicroTriangle::operator==(const MicroTriangle& triangle)
{
    uint32_t i, j;
    bool isFinded;

    for(i=0; i<3; i++) if(_nodes[i] == EmptyNode()) return false;
    for(i=0; i<3; i++) if(triangle._nodes[i] == EmptyNode()) return false;
    for(i=0; i<3; i++)
    {
        isFinded = false;
        for(j=0; j<3; j++)
            if(_nodes[i] == triangle._nodes[j])
            {
                isFinded = true;
                break;
            }
        if(!isFinded) return false;
    }
    return true;
}

/***********************************************************************************
 *                    Другие функции                                               *
 ***********************************************************************************/

bool IsSharedTetraedrs(MicroTetraedr& tetraedr1, MicroTetraedr& tetraedr2)
{
    uint32_t i, count = 0;
    for(i=0; i<4; i++) if(tetraedr1.isNodeTetraedr(tetraedr2.nodes(i)->point())) count++;
    if(count == 3) return true;
    else return false;
}

bool IsCrossingTetraedrs(MicroTetraedr &tetraedr1, MicroTetraedr &tetraedr2)
{
    uint32_t i;
    uint32_t count;

    count = 0;
    for(i=0; i<4; i++)
    {
        if(tetraedr1.isInTetraedr(tetraedr2.nodes(i)->point()))
        {
            if(tetraedr1.isNodeTetraedr(tetraedr2.nodes(i)->point()))
                count++;
            else return true;
        }
    }
    if(count == 4) return true;
    count = 0;
    for(i=0; i<4; i++)
        if(tetraedr2.isInTetraedr(tetraedr1.nodes(i)->point()))
            if(!tetraedr2.isNodeTetraedr(tetraedr1.nodes(i)->point())) return true;

    return false;
}

bool IsDeloneTetraedrs(ListMicroTetraedr::iterator &tetraedr1, ListMicroTetraedr::iterator &tetraedr2)
{
    uint32_t j;
    if(tetraedr1 == EmptyIterator() || tetraedr2 == EmptyIterator()) return true;
    if(!IsSharedTetraedrs(*tetraedr1, *tetraedr2)) return false;

    for(j=0; j<4; j++) if(tetraedr1->tetraedrs(j) == tetraedr2) break;   /* поиск указателя на себя */
#ifdef DEBUG_MODE
    if(j == 4)
    {
        printf("Error in find pointer in function IsDeloneTetraedrs \n");
        exit(-1);
    }
#endif
    if(tetraedr1->isInSphere((tetraedr2->nodes(j + 3))->point())) return false;

    for(j=0; j<4; j++) if(tetraedr2->tetraedrs(j) == tetraedr1) break;   /* поиск указателя на себя */
#ifdef DEBUG_MODE
    if(j == 4)
    {
        printf("Error in find pointer in function IsDeloneTetraedrs \n");
        exit(-1);
    }
#endif
    if(tetraedr2->isInSphere((tetraedr1->nodes(j + 3))->point())) return false;

    return true;
}

bool IsDeloneTetraedr(ListMicroTetraedr::iterator &tetraedr)
{
    uint32_t i;
    for(i=0; i<4; i++)
        if(!IsDeloneTetraedrs(tetraedr, tetraedr->tetraedrs(i))) return false;
    return true;
}

bool IsSingularTetraedr(MathPoint3D &point1, MathPoint3D &point2, MathPoint3D &point3, MathPoint3D &point4, double epsilonVolume)
{
    if(VolumeTetraedr(point1, point2, point3, point4) < epsilonVolume) return true;
    else return false;
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


//EOF
