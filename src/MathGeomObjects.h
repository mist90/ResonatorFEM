#ifndef MATHGEOMOBJECTS
#define MATHGEOMOBJECTS

#include "MathPoint3D.h"
#include "MathVector3D.h"

MathPoint3D GetPointFromCircle(double rad, const MathVector3D normal, double angleXY);

MathPoint3D IntersectionLinesSurface(const MathPoint3D point1Lines,
                                     const MathPoint3D point2Lines,
                                     const MathPoint3D point1Surface,
                                     const MathPoint3D point2Surface,
                                     const MathPoint3D point3Surface);

MathPoint3D ProjectionPointSurface(const MathPoint3D point,
                                   const MathPoint3D point1Surface,
                                   const MathPoint3D point2Surface,
                                   const MathPoint3D point3Surface);

/* площадь треугольника */
double AreaTriangle(const MathPoint3D& point1, const MathPoint3D& point2, const MathPoint3D& point3);

#endif // MATHGEOMOBJECTS

