#include "MathVector3D.h"

MathVector3D::MathVector3D(const double &x, const double &y, const double &z)
{
    _x = x;
    _y = y;
    _z = z;
    _beginPoint = MathPoint3D(0.0, 0.0, 0.0);
}

MathVector3D::MathVector3D(const double &x, const double &y, const double &z, const MathPoint3D &begin)
{
    _x = x;
    _y = y;
    _z = z;
    _beginPoint = begin;
}

MathVector3D::MathVector3D(const MathPoint3D &begin, const MathPoint3D &end)
{
    _x = ((MathPoint3D)end).getX() - ((MathPoint3D)begin).getX();
    _y = ((MathPoint3D)end).getY() - ((MathPoint3D)begin).getY();
    _z = ((MathPoint3D)end).getZ() - ((MathPoint3D)begin).getZ();
    _beginPoint = begin;
}

MathVector3D::MathVector3D()
{
    _x = _y = _z = 0;
    _beginPoint = MathPoint3D(0.0, 0.0, 0.0);
}

MathVector3D::MathVector3D(const MathVector3D& p)
{
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _beginPoint = p._beginPoint;
}

void MathVector3D::operator=(const MathVector3D& p)
{
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _beginPoint = p._beginPoint;
}

MathVector3D MathVector3D::operator+(const MathVector3D& p)
{
    MathVector3D s;
    s._x = _x + p._x;
    s._y = _y + p._y;
    s._z = _z + p._z;
    s._beginPoint = _beginPoint;
    return s;
}

MathVector3D MathVector3D::operator-(const MathVector3D& p)
{
    MathVector3D s;
    s._x = _x - p._x;
    s._y = _y - p._y;
    s._z = _z - p._z;
    s._beginPoint = _beginPoint;
    return s;
}

MathVector3D MathVector3D::operator-()
{
    MathVector3D s;
    s._x = -_x;
    s._y = -_y;
    s._z = -_z;
    s._beginPoint = _beginPoint;
    return s;
}

double MathVector3D::operator*(const MathVector3D& p)
{
    return _x * p._x + _y * p._y + _z * p._z;
}

MathVector3D MathVector3D::operator*(const double& mul)
{
    MathVector3D ret;
    ret.setX(getX()*mul);
    ret.setY(getY()*mul);
    ret.setZ(getZ()*mul);
    ret.setBegin(getBegin());
    return ret;
}

MathVector3D MathVector3D::operator ^(const MathVector3D& p)
{
    MathVector3D ret;

    ret.setX(getY()*((MathVector3D)p).getZ() - getZ()*((MathVector3D)p).getY());
    ret.setY(getZ()*((MathVector3D)p).getX() - getX()*((MathVector3D)p).getZ());
    ret.setZ(getX()*((MathVector3D)p).getY() - getY()*((MathVector3D)p).getX());
    ret.setBegin(getBegin());
    return ret;
}

void MathVector3D::setX(const double &x)
{
    _x = x;
}

void MathVector3D::setY(const double &y)
{
    _y = y;
}

void MathVector3D::setZ(const double &z)
{
    _z = z;
}

void MathVector3D::setBegin(const MathPoint3D &point)
{
    _beginPoint = point;
}

double MathVector3D::getX()
{
    return _x;
}

double MathVector3D::getY()
{
    return _y;
}

double MathVector3D::getZ()
{
    return _z;
}

MathPoint3D MathVector3D::getBegin()
{
    return _beginPoint;
}

MathPoint3D MathVector3D::getEnd()
{
    return MathPoint3D(_beginPoint.getX() + _x, _beginPoint.getY() + _y, _beginPoint.getZ() + _z);
}

double MathVector3D::lenght()
{
    return sqrt(_x*_x + _y*_y + _z*_z);
}

MathVector3D MathVector3D::unitVector()
{
    MathVector3D unit;
    double len = lenght();

    if(!len) len = 1.0;
    unit._x = _x/len;
    unit._y = _y/len;
    unit._z = _z/len;
    unit._beginPoint = _beginPoint;
    return unit;
}

//EOF

