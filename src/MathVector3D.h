#ifndef MATHVECTOR3D_H
#define MATHVECTOR3D_H
#include <math.h>
#include "MathPoint3D.h"

class MathVector3D
{
public:
    MathVector3D(const double& x, const double& y, const double& z);
    MathVector3D(const double& x, const double& y, const double& z, const MathPoint3D& begin);
    MathVector3D(const MathPoint3D& begin, const MathPoint3D& end);
    MathVector3D();
    MathVector3D(const MathVector3D& p);
    void         operator=(const MathVector3D& p);
    MathVector3D operator+(const MathVector3D& p);
    MathVector3D operator-(const MathVector3D& p);
    MathVector3D operator-();
    double       operator*(const MathVector3D& p);  // скалярное произведение
    MathVector3D operator*(const double& mul);
    MathVector3D operator^(const MathVector3D& p);  // векторное произведение
    void         setX(const double& x);
    void         setY(const double& y);
    void         setZ(const double& z);
    void         setBegin(const MathPoint3D& point);
    double       getX();
    double       getY();
    double       getZ();
    MathPoint3D  getBegin();
    MathPoint3D  getEnd();
    double       lenght();
    MathVector3D unitVector();
private:
    MathPoint3D  _beginPoint;
    double       _x;
    double       _y;
    double       _z;
};


#endif // MATHVECTOR3D_H
