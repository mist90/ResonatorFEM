#ifndef MATHOBJECT
#define MATHOBJECT
#include "MathVector3D.h"
#include "MathPoint3D.h"
#include "MathGeomObjects.h"
#include <vector>
#include <stdint.h>

#define EPSILON_MUL (0.1)

#define AUTO_STEP (-1)

#define AXIS_X  (0)
#define AXIS_Y  (1)
#define AXIS_Z  (2)

#define FIGURE_NOT_DEFINED     (0)
#define FIGURE_SPHERE          (1)
#define FIGURE_CILINDR         (2)
#define FIGURE_PARALLELEPIPED  (3)
#define FIGURE_SPIRAL          (4)


class MathObject
{
public:
    MathObject();
    void makeSphere(const MathPoint3D& center,
                    const double& rad,
                    const double& step = AUTO_STEP);
    void makeCilindr(const MathPoint3D& center,
                     const double& rad,
                     const double& len,
                     const double& step = AUTO_STEP);
    void makeParallelepiped(const MathPoint3D& center,
                     const double& deltaX,
                     const double& deltaY,
                     const double& deltaZ,
                     const double& step = AUTO_STEP);
    void makeSpiral(const MathPoint3D& center,
                     const double& radDepth,
                     const double& radSpiral,
                     const double& stepSpiral,
                     const double& angleBegin,
                     const double &lenSpiral,
                     const double& step = AUTO_STEP);
    bool rotate(uint32_t axis,const double& angle);
    void setVector(const MathVector3D& normal);
    bool getPoints(std::vector<MathPoint3D>& points, bool isSolid);
    bool isPointInFigure(const MathPoint3D& point);
    MathPoint3D translateBackPoint(const MathPoint3D& point);
    void clear();

private:
    bool getPointsCilindr(const MathPoint3D& center,
                          const double& rad,
                          const double& len,
                          const double& step,
                          std::vector<MathPoint3D>& points,
                          bool isSolid);
    bool getPointsSphere(const MathPoint3D& center,
                         const double& rad,
                         const double& step,
                         std::vector<MathPoint3D>& points,
                         bool isSolid);
    bool getPointsSpiral(const MathPoint3D& center,
                         const double& radDepth,
                         const double& radSpiral,
                         const double& stepSpiral,
                         const double& angleBegin,
                         const double &lenSpiral,
                         const double& step,
                         std::vector<MathPoint3D>& points,
                         bool isSolid);
    bool isPointInSpiral(const MathPoint3D &point);
    uint32_t typeFigure;
    std::vector<double> vectorParameters;
    std::vector<double> vectorValueRotates;
    std::vector<uint32_t> vectorTypeRotates;
    double sinNormal[2], cosNormal[2];
    MathPoint3D centerFigure;
};

#endif // MATHOBJECT

