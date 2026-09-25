#ifndef MATHVECTOR3D_H
#define MATHVECTOR3D_H
#include <math.h>

/* A 3D vector, also used as a point (position vector). */
class MathVector3D
{
public:
    MathVector3D(const double& x, const double& y, const double& z);
    MathVector3D();
    MathVector3D(const MathVector3D& p);
    void         operator=(const MathVector3D& p);
    bool         operator==(const MathVector3D& p);
    MathVector3D operator+(const MathVector3D& p);
    MathVector3D operator-(const MathVector3D& p);
    MathVector3D operator-();
    double       operator*(const MathVector3D& p);  /* dot product   */
    MathVector3D operator*(const double& mul);
    MathVector3D operator/(const double& value);
    MathVector3D operator^(const MathVector3D& p);  /* cross product */
    void         setX(const double& x);
    void         setY(const double& y);
    void         setZ(const double& z);
    double       getX();
    double       getY();
    double       getZ();
    double       lenght();
    MathVector3D unitVector();
private:
    double       _x;
    double       _y;
    double       _z;
};

#endif // MATHVECTOR3D_H
