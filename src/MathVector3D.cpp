#include "MathVector3D.h"

MathVector3D::MathVector3D(const double &x, const double &y, const double &z)
{
    _x = x;
    _y = y;
    _z = z;
}

MathVector3D::MathVector3D()
{
    _x = _y = _z = 0;
}

MathVector3D::MathVector3D(const MathVector3D& p)
{
    _x = p._x;
    _y = p._y;
    _z = p._z;
}

void MathVector3D::operator=(const MathVector3D& p)
{
    _x = p._x;
    _y = p._y;
    _z = p._z;
}

bool MathVector3D::operator==(const MathVector3D& p)
{
    return _x == p._x && _y == p._y && _z == p._z;
}

MathVector3D MathVector3D::operator+(const MathVector3D& p)
{
    return MathVector3D(_x + p._x, _y + p._y, _z + p._z);
}

MathVector3D MathVector3D::operator-(const MathVector3D& p)
{
    return MathVector3D(_x - p._x, _y - p._y, _z - p._z);
}

MathVector3D MathVector3D::operator-()
{
    return MathVector3D(-_x, -_y, -_z);
}

double MathVector3D::operator*(const MathVector3D& p)
{
    return _x * p._x + _y * p._y + _z * p._z;
}

MathVector3D MathVector3D::operator*(const double& mul)
{
    return MathVector3D(_x * mul, _y * mul, _z * mul);
}

MathVector3D MathVector3D::operator/(const double& value)
{
    return MathVector3D(_x / value, _y / value, _z / value);
}

MathVector3D MathVector3D::operator ^(const MathVector3D& p)
{
    return MathVector3D(_y * p._z - _z * p._y,
                        _z * p._x - _x * p._z,
                        _x * p._y - _y * p._x);
}

void MathVector3D::setX(const double &x) { _x = x; }
void MathVector3D::setY(const double &y) { _y = y; }
void MathVector3D::setZ(const double &z) { _z = z; }

double MathVector3D::getX() { return _x; }
double MathVector3D::getY() { return _y; }
double MathVector3D::getZ() { return _z; }

double MathVector3D::lenght()
{
    return sqrt(_x*_x + _y*_y + _z*_z);
}

MathVector3D MathVector3D::unitVector()
{
    double len = lenght();
    if(!len) len = 1.0;
    return MathVector3D(_x/len, _y/len, _z/len);
}
