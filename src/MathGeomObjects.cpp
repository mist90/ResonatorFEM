#include "MathGeomObjects.h"

MathPoint3D GetPointFromCircle(double rad, const MathVector3D normal, double angleXY)
{
    MathPoint3D ret(rad, 0.0, 0.0);
    MathVector3D vector = normal;
    double divXY;

    vector = vector.unitVector();
    divXY = sqrt(vector.getX()*vector.getX() + vector.getY()*vector.getY());
    ret = ret.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), angleXY);
    ret = ret.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), sqrt(1 - vector.getZ()*vector.getZ()), vector.getZ());
    if(divXY != 0.0) ret = ret.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), vector.getY()/divXY, vector.getX()/divXY);
    ret = ret + vector.getBegin();
    return ret;
}

MathPoint3D IntersectionLinesSurface(const MathPoint3D point1Lines,const MathPoint3D point2Lines,const MathPoint3D point1Surface,const MathPoint3D point2Surface,const MathPoint3D point3Surface)
{
    MathPoint3D ret;
    MathVector3D vec1(point1Surface, point2Surface), vec2(point1Surface, point3Surface), vecNormal, vecLine;
    double dt;

    vecNormal = vec1^vec2;
    vecLine = MathVector3D(point1Lines, point2Lines);
    vec1 = MathVector3D(point1Surface, point1Lines);
    vec2 = MathVector3D(point1Surface, point2Lines);
    dt = -(vecNormal*vec1)/(vecNormal*vecLine);
    ret = (((MathPoint3D)point1Lines) + MathPoint3D(vecLine.getX()*dt, vecLine.getY()*dt, vecLine.getZ()*dt));
    return ret;
}

MathPoint3D ProjectionPointSurface(const MathPoint3D point, const MathPoint3D point1Surface, const MathPoint3D point2Surface, const MathPoint3D point3Surface)
{
    MathPoint3D ret;
    MathVector3D vec1(point1Surface, point2Surface), vec2(point1Surface, point3Surface), vecNormal;
    double dt;

    vecNormal = vec1^vec2;
    vec1 = MathVector3D(point1Surface, point);
    dt = -(vecNormal*vec1)/(vecNormal*vecNormal);
    ret = ((MathPoint3D)point) + MathPoint3D(vecNormal.getX()*dt, vecNormal.getY()*dt, vecNormal.getZ()*dt);
    return ret;
}


double AreaTriangle(const MathPoint3D &point1, const MathPoint3D &point2, const MathPoint3D &point3)
{
    double sx, sy, sz;
    MathPoint3D *_point1 = (MathPoint3D*)&point1, *_point2 = (MathPoint3D*)&point2, *_point3 = (MathPoint3D*)&point3;

    sx = (_point2->getY() - _point1->getY())*(_point3->getZ() - _point1->getZ()) - (_point3->getY() - _point1->getY())*(_point2->getZ() - _point1->getZ());
    sy = -((_point2->getX() - _point1->getX())*(_point3->getZ() - _point1->getZ()) - (_point3->getX() - _point1->getX())*(_point2->getZ() - _point1->getZ()));
    sz = (_point2->getX() - _point1->getX())*(_point3->getY() - _point1->getY()) - (_point3->getX() - _point1->getX())*(_point2->getY() - _point1->getY());
    return sqrt(sx*sx + sy*sy +sz*sz)/2.0;
}


//EOF
