#ifndef MATHPOINT3D_H
#define MATHPOINT3D_H


/* класс точка с тремя коордтнатами */
class MathPoint3D
{
public:
    MathPoint3D(const double& x, const double& y, const double& z);
    MathPoint3D();
    MathPoint3D(const MathPoint3D& p);
    void operator=(const MathPoint3D& p);
    bool operator==(const MathPoint3D& p);
    MathPoint3D operator+(const MathPoint3D& p);
    MathPoint3D operator-(const MathPoint3D& p);
    MathPoint3D operator*(const double& value);
    MathPoint3D operator/(const double& value);
    void setX(const double& x);
    void setY(const double& y);
    void setZ(const double& z);
    double getX();
    double getY();
    double getZ();
    MathPoint3D  rotateParallelAxisX(const MathPoint3D& pointCenter, const double& angle);
    MathPoint3D  rotateParallelAxisX(const MathPoint3D& pointCenter, const double& sinAngle, const double& cosAngle);
    MathPoint3D  rotateParallelAxisY(const MathPoint3D& pointCenter, const double& angle);
    MathPoint3D  rotateParallelAxisY(const MathPoint3D& pointCenter, const double& sinAngle, const double& cosAngle);
    MathPoint3D  rotateParallelAxisZ(const MathPoint3D& pointCenter, const double& angle);
    MathPoint3D  rotateParallelAxisZ(const MathPoint3D& pointCenter, const double& sinAngle, const double& cosAngle);
    friend double LenSquaredPoints(const MathPoint3D& point1, const MathPoint3D& point2);
    friend double LenPoints(const MathPoint3D& point1, const MathPoint3D& point2);
private:
    double _x;
    double _y;
    double _z;
};

/* Возвращает расстояние в квадрате между точками */
double LenSquaredPoints(const MathPoint3D& point1, const MathPoint3D& point2);

/* Возвращает расстояние между точками */
double LenPoints(const MathPoint3D& point1, const MathPoint3D& point2);


#endif // MATHPOINT3D_H
