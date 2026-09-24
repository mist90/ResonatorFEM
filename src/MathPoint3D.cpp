#include "MathPoint3D.h"
#include <math.h>

MathPoint3D::MathPoint3D(const double &x, const double &y, const double &z)
{
    _x = x;
    _y = y;
    _z = z;

}

MathPoint3D::MathPoint3D()
{
    _x = _y = _z = 0;
}

MathPoint3D::MathPoint3D(const MathPoint3D& p)
{
    _x = p._x;
    _y = p._y;
    _z = p._z;
}

void MathPoint3D::operator=(const MathPoint3D& p)
{
    _x = p._x;
    _y = p._y;
    _z = p._z;
}

bool MathPoint3D::operator==(const MathPoint3D& p)
{
    if(_x != p._x) return false;
    if(_y != p._y) return false;
    if(_z != p._z) return false;
    return true;
}

MathPoint3D MathPoint3D::operator +(const MathPoint3D& p)
{
    MathPoint3D ret;
    ret.setX(_x + p._x);
    ret.setY(_y + p._y);
    ret.setZ(_z + p._z);
    return ret;
}

MathPoint3D MathPoint3D::operator -(const MathPoint3D& p)
{
    MathPoint3D ret;
    ret.setX(_x - p._x);
    ret.setY(_y - p._y);
    ret.setZ(_z - p._z);
    return ret;
}

MathPoint3D MathPoint3D::operator *(const double& value)
{
    MathPoint3D ret;
    ret.setX(_x * value);
    ret.setY(_y * value);
    ret.setZ(_z * value);
    return ret;
}

MathPoint3D MathPoint3D::operator /(const double& value)
{
    MathPoint3D ret;
    ret.setX(_x / value);
    ret.setY(_y / value);
    ret.setZ(_z / value);
    return ret;
}

void MathPoint3D::setX(const double &x)
{
    _x = x;
}

void MathPoint3D::setY(const double &y)
{
    _y = y;
}

void MathPoint3D::setZ(const double &z)
{
    _z = z;
}

double MathPoint3D::getX()
{
    return _x;
}

double MathPoint3D::getY()
{
    return _y;
}

double MathPoint3D::getZ()
{
    return _z;
}

MathPoint3D MathPoint3D::rotateParallelAxisX(const MathPoint3D &pointCenter, const double &angle)
{
    double angleRad = angle*M_PI/180.0;
    MathPoint3D ret;

    ret.setX(_x);
    ret.setY(((MathPoint3D)pointCenter).getY() + (_y - ((MathPoint3D)pointCenter).getY())*cos(angleRad) - (_z - ((MathPoint3D)pointCenter).getZ())*sin(angleRad));
    ret.setZ(((MathPoint3D)pointCenter).getZ() + (_y - ((MathPoint3D)pointCenter).getY())*sin(angleRad) + (_z - ((MathPoint3D)pointCenter).getZ())*cos(angleRad));
    return ret;
}

MathPoint3D MathPoint3D::rotateParallelAxisX(const MathPoint3D &pointCenter, const double &sinAngle, const double &cosAngle)
{
    MathPoint3D ret;

    ret.setX(_x);
    ret.setY(((MathPoint3D)pointCenter).getY() + (_y - ((MathPoint3D)pointCenter).getY())*cosAngle - (_z - ((MathPoint3D)pointCenter).getZ())*sinAngle);
    ret.setZ(((MathPoint3D)pointCenter).getZ() + (_y - ((MathPoint3D)pointCenter).getY())*sinAngle + (_z - ((MathPoint3D)pointCenter).getZ())*cosAngle);
    return ret;
}

MathPoint3D MathPoint3D::rotateParallelAxisY(const MathPoint3D &pointCenter, const double &angle)
{
    double angleRad = angle*M_PI/180.0;
    MathPoint3D ret;

    ret.setX(((MathPoint3D)pointCenter).getX() + (_z - ((MathPoint3D)pointCenter).getZ())*sin(angleRad) + (_x - ((MathPoint3D)pointCenter).getX())*cos(angleRad));
    ret.setY(_y);
    ret.setZ(((MathPoint3D)pointCenter).getZ() + (_z - ((MathPoint3D)pointCenter).getZ())*cos(angleRad) - (_x - ((MathPoint3D)pointCenter).getX())*sin(angleRad));
    return ret;
}

MathPoint3D MathPoint3D::rotateParallelAxisY(const MathPoint3D &pointCenter, const double &sinAngle, const double &cosAngle)
{
    MathPoint3D ret;

    ret.setX(((MathPoint3D)pointCenter).getX() + (_z - ((MathPoint3D)pointCenter).getZ())*sinAngle + (_x - ((MathPoint3D)pointCenter).getX())*cosAngle);
    ret.setY(_y);
    ret.setZ(((MathPoint3D)pointCenter).getZ() + (_z - ((MathPoint3D)pointCenter).getZ())*cosAngle - (_x - ((MathPoint3D)pointCenter).getX())*sinAngle);
    return ret;
}

MathPoint3D MathPoint3D::rotateParallelAxisZ(const MathPoint3D &pointCenter, const double &angle)
{
    double angleRad = angle*M_PI/180.0;
    MathPoint3D ret;

    ret.setX(((MathPoint3D)pointCenter).getX() + (_x - ((MathPoint3D)pointCenter).getX())*cos(angleRad) - (_y - ((MathPoint3D)pointCenter).getY())*sin(angleRad));
    ret.setY(((MathPoint3D)pointCenter).getY() + (_x - ((MathPoint3D)pointCenter).getX())*sin(angleRad) + (_y - ((MathPoint3D)pointCenter).getY())*cos(angleRad));
    ret.setZ(_z);
    return ret;
}

MathPoint3D MathPoint3D::rotateParallelAxisZ(const MathPoint3D &pointCenter, const double &sinAngle, const double &cosAngle)
{
    MathPoint3D ret;

    ret.setX(((MathPoint3D)pointCenter).getX() + (_x - ((MathPoint3D)pointCenter).getX())*cosAngle - (_y - ((MathPoint3D)pointCenter).getY())*sinAngle);
    ret.setY(((MathPoint3D)pointCenter).getY() + (_x - ((MathPoint3D)pointCenter).getX())*sinAngle + (_y - ((MathPoint3D)pointCenter).getY())*cosAngle);
    ret.setZ(_z);
    return ret;
}

double LenSquaredPoints(const MathPoint3D& point1, const MathPoint3D& point2)
{
    return (point2._x - point1._x)*(point2._x - point1._x) + (point2._y - point1._y)*(point2._y - point1._y) + (point2._z - point1._z)*(point2._z - point1._z);
}

double LenPoints(const MathPoint3D& point1, const MathPoint3D& point2)
{
    return sqrt((point2._x - point1._x)*(point2._x - point1._x) + (point2._y - point1._y)*(point2._y - point1._y) + (point2._z - point1._z)*(point2._z - point1._z));
}


//EOF
