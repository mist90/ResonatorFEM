#include "MathObject.h"

MathObject::MathObject()
{
    typeFigure = FIGURE_NOT_DEFINED;
    centerFigure = MathPoint3D(0.0, 0.0, 0.0);
}

void MathObject::makeSphere(const MathPoint3D &center, const double &rad, const double &step)
{
    double _step;

    clear();
    if(step != AUTO_STEP && step < 0) return;
    if(step == AUTO_STEP)
        _step = 2*M_PI*rad/24.0;
    else _step = 2*M_PI*rad/24.0*step;
    typeFigure = FIGURE_SPHERE;
    centerFigure = center;
    vectorParameters.push_back(fabs(rad));
    vectorParameters.push_back(_step);
}

void MathObject::makeCilindr(const MathPoint3D &center,
                             const double &rad,
                             const double &len,
                             const double &step)
{
    double _step;

    clear();
    if(step != AUTO_STEP && step < 0) return;
    if(step == AUTO_STEP)
        _step = 2*M_PI*rad/24.0;
    else _step = 2*M_PI*rad/24.0*step;
    typeFigure = FIGURE_CILINDR;
    centerFigure = center;
    vectorParameters.push_back(fabs(rad));
    vectorParameters.push_back(len);
    vectorParameters.push_back(_step);
}

void MathObject::makeParallelepiped(const MathPoint3D &center,
                                    const double &deltaX,
                                    const double &deltaY,
                                    const double &deltaZ,
                                    const double &step)
{
    clear();
    if(step != AUTO_STEP && step < 0) return;

    typeFigure = FIGURE_PARALLELEPIPED;
    centerFigure = center;
    vectorParameters.push_back(fabs(deltaX));
    vectorParameters.push_back(fabs(deltaY));
    vectorParameters.push_back(fabs(deltaZ));
    vectorParameters.push_back(step);
}

void MathObject::makeSpiral(const MathPoint3D &center,
                            const double &radDepth,
                            const double &radSpiral,
                            const double &stepSpiral,
                            const double &angleBegin,
                            const double &lenSpiral,
                            const double &step)
{
    double _step;

    clear();
    if(step != AUTO_STEP && step < 0) return;
    if(step == AUTO_STEP)
        _step = 2*M_PI*radDepth/12.0;
    else _step = 2*M_PI*radDepth/12.0*step;
    typeFigure = FIGURE_SPIRAL;
    centerFigure = center;
    vectorParameters.push_back(fabs(radDepth));
    vectorParameters.push_back(fabs(radSpiral));
    vectorParameters.push_back(fabs(stepSpiral));
    vectorParameters.push_back(angleBegin);
    vectorParameters.push_back(fabs(lenSpiral));
    vectorParameters.push_back(_step);
}

bool MathObject::rotate(uint32_t axis, const double &angle)
{
    if(axis == AXIS_X || axis == AXIS_Y || axis == AXIS_Z)
    {
        vectorTypeRotates.push_back(axis);
        vectorValueRotates.push_back(angle);
        return true;
    }
    else return false;
}

void MathObject::setVector(const MathVector3D &normal)
{
    vectorTypeRotates.clear();
    vectorValueRotates.clear();
    MathVector3D vector = normal;
    double divXY;

    vector = vector.unitVector();
    divXY = sqrt(vector.getX()*vector.getX() + vector.getY()*vector.getY());
    sinNormal[0] = sqrt(1 - vector.getZ()*vector.getZ());
    cosNormal[0] = vector.getZ();
    if(divXY == 0.0)
    {
        sinNormal[1] = 0.0;
        cosNormal[1] = 1.0;
    }
    else
    {
        sinNormal[1] = vector.getY()/divXY;
        cosNormal[1] = vector.getX()/divXY;
    }
}

bool MathObject::getPoints(std::vector<MathPoint3D> &points, bool isSolid)
{
    double rad = vectorParameters[0];
    double radStep, currentStep;
    uint32_t i, numCilindr;

    switch(typeFigure)
    {
        case FIGURE_CILINDR:
        if(!isSolid)
            return getPointsCilindr(centerFigure,
                    vectorParameters[0],
                    vectorParameters[1],
                    vectorParameters[2]/2.0,
                    points,
                    isSolid);
        else
        {
            radStep = vectorParameters[2];
            numCilindr = rad/radStep;
            radStep = rad/numCilindr;
            for(i=0; i<numCilindr; i++)
            {
                if(i == numCilindr - 1) currentStep = vectorParameters[2]/2.0;
                else currentStep = vectorParameters[2];
                if(!getPointsCilindr(centerFigure,
                                     (i + 1)*radStep,
                                     vectorParameters[1],
                                     currentStep,
                                     points,
                                     isSolid)) return false;
            }
            return true;
        }
            break;

        case FIGURE_SPHERE:
            return getPointsSphere(centerFigure,
                    vectorParameters[0],
                    vectorParameters[1],
                    points,
                    isSolid);
            break;

        case FIGURE_SPIRAL:
            return getPointsSpiral(centerFigure,
                    vectorParameters[0],
                    vectorParameters[1],
                    vectorParameters[2],
                    vectorParameters[3],
                    vectorParameters[4],
                    vectorParameters[5],
                    points,
                    isSolid);
            break;
    }

    return true;
}

void MathObject::clear()
{
    typeFigure = FIGURE_NOT_DEFINED;
    centerFigure = MathPoint3D(0.0, 0.0, 0.0);
    vectorParameters.clear();
    vectorTypeRotates.clear();
    vectorValueRotates.clear();
    sinNormal[0] = sinNormal[1] = 0.0;
    cosNormal[0] = cosNormal[1] = 1.0;
}

bool MathObject::getPointsCilindr(const MathPoint3D &center,
                                  const double &rad,
                                  const double &len,
                                  const double &step,
                                  std::vector<MathPoint3D> &points,
                                  bool isSolid)
{
    uint32_t i, j, k;
    uint32_t numLen, numRad;
    double deltaLen, deltaRad;
    MathPoint3D point;

    /* Настройка шагов */

    numLen = len/step;
    if(!numLen) return false;
    deltaLen = len/numLen;
    deltaRad = 4*step/rad;
    numRad = 2*M_PI/deltaRad;
    if(!numRad) return false;
    deltaRad = 2*M_PI/numRad;

    /* Генерация точек оси */
    if(isSolid)
        for(i=0; i<=numLen/2; i++)
        {
            point.setX(0.0);
            point.setY(0.0);
            point.setZ(2*deltaLen*i - len/2.0);
            point = point.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), sinNormal[0], cosNormal[0]);
            point = point.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), sinNormal[1], cosNormal[1]);
            for(k=0; k<vectorTypeRotates.size(); k++)
                switch(vectorTypeRotates[k])
                {
                case AXIS_X:
                    point = point.rotateParallelAxisX(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;

                case AXIS_Y:
                    point = point.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;

                case AXIS_Z:
                    point = point.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;
                }

            point = point + center;
            points.push_back(point);
        }

    for(i=0; i<=numLen; i++)
    {
        for(j=0; j<numRad; j++)
        {
            point.setX(rad*cos(j*deltaRad + i*deltaRad/2.0));
            point.setY(rad*sin(j*deltaRad + i*deltaRad/2.0));
            point.setZ(deltaLen*i - len/2.0);
            point = point.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), sinNormal[0], cosNormal[0]);
            point = point.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), sinNormal[1], cosNormal[1]);
            for(k=0; k<vectorTypeRotates.size(); k++)
                switch(vectorTypeRotates[k])
                {
                case AXIS_X:
                    point = point.rotateParallelAxisX(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;

                case AXIS_Y:
                    point = point.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;

                case AXIS_Z:
                    point = point.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;
                }
            point = point + center;
            points.push_back(point);
        }
    }
    return true;
}


bool MathObject::getPointsSphere(const MathPoint3D &center, const double &rad, const double &step, std::vector<MathPoint3D> &points, bool isSolid)
{
    uint32_t i, j;
    uint32_t numRad;
    double deltaRad;
    MathPoint3D point;


    deltaRad = step/rad;
    numRad = 2*M_PI/deltaRad;
    if(!numRad) return false;
    deltaRad = 2*M_PI/numRad;

    if(isSolid) points.push_back(center);

    for(i=0; i<numRad; i++)
    {
        for(j=0; j<numRad; j++)
        {
            point.setX(rad*cos(j*deltaRad)*cos(i*deltaRad));
            point.setY(rad*sin(j*deltaRad));
            point.setZ(rad*cos(j*deltaRad)*sin(i*deltaRad));
            point = point + center;
            points.push_back(point);
        }
    }
    return true;
}

bool MathObject::getPointsSpiral(const MathPoint3D &center,
                                 const double &radDepth,
                                 const double &radSpiral,
                                 const double &stepSpiral,
                                 const double &angleBegin,
                                 const double &lenSpiral,
                                 const double &step,
                                 std::vector<MathPoint3D> &points,
                                 bool isSolid)
{
    uint32_t j, k;
    uint32_t numRad;
    double deltaRad;
    double numStep;
    MathPoint3D point, point2;
    MathVector3D n;

    deltaRad = 4*step/radDepth;
    numRad = 2*M_PI/deltaRad;
    if(numRad < 4)
    {
        numRad = 4;
        deltaRad = 2*M_PI/numRad;
    }
    deltaRad = 2*M_PI/numRad;

    numStep = 0;
    while(numStep*stepSpiral/(2*M_PI*radSpiral) <= lenSpiral)
    {
        point.setX(radSpiral*cos(numStep/(radSpiral) + angleBegin*M_PI/180.0));
        point.setY(radSpiral*sin(numStep/(radSpiral) + angleBegin*M_PI/180.0));
        point.setZ(numStep*stepSpiral/(2*M_PI*radSpiral) - lenSpiral/2.0);
        n = MathVector3D(-sin(numStep/(radSpiral) + angleBegin*M_PI/180.0), cos(numStep/(radSpiral) + angleBegin*M_PI/180.0), stepSpiral/(2*M_PI*radSpiral), point);

        for(j=0; j<numRad; j++)
        {
            point2 = GetPointFromCircle(radDepth, n, (j*deltaRad + numStep*deltaRad/step/2.0)*180.0/M_PI);
            point2 = point2.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), sinNormal[0], cosNormal[0]);
            point2 = point2.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), sinNormal[1], cosNormal[1]);
            for(k=0; k<vectorTypeRotates.size(); k++)
                switch(vectorTypeRotates[k])
                {
                case AXIS_X:
                    point2 = point2.rotateParallelAxisX(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;

                case AXIS_Y:
                    point2 = point2.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;

                case AXIS_Z:
                    point2 = point2.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                    break;
                }
            point2 = point2 + center;
            points.push_back(point2);
        }
        point = point.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), sinNormal[0], cosNormal[0]);
        point = point.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), sinNormal[1], cosNormal[1]);
        for(k=0; k<vectorTypeRotates.size(); k++)
            switch(vectorTypeRotates[k])
            {
            case AXIS_X:
                point = point.rotateParallelAxisX(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                break;

            case AXIS_Y:
                point = point.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                break;

            case AXIS_Z:
                point = point.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), vectorValueRotates[k]);
                break;
            }
        point = point + center;
        if(isSolid) points.push_back(point);
        numStep += step;
    }
    return true;
}

bool MathObject::isPointInSpiral(const MathPoint3D &point)
{
    double angle;
    double rad;
    double radDepth = vectorParameters[0];
    double radSpiral = vectorParameters[1];
    double stepSpiral = vectorParameters[2];
    double angleBegin = vectorParameters[3]*M_PI/180.0;
    double len = vectorParameters[4];
    MathPoint3D pointCenterCircle, _point = point;
    double zBegin, zPeriod;
    int32_t nPeriod;

    /* Нахождение центра окружности, перпендикулярной спирали и содержащей точку point */
    if(_point.getX() != 0.0 || _point.getY() != 0.0)
        angle = asin(_point.getY()/sqrt(_point.getX()*_point.getX() + _point.getY()*_point.getY()));
    else angle = 0.0;
    if(_point.getX() < 0.0)
    {
        if(angle > 0.0) angle = M_PI - angle;
        else angle = -(M_PI + angle);
    }
    pointCenterCircle.setX(radSpiral*cos(angle));
    pointCenterCircle.setY(radSpiral*sin(angle));
    /* M_PI/stepSpiral*len появляется из-за смещения центра отсчета в середину спирали */
    zBegin = stepSpiral*(angle - angleBegin - M_PI/stepSpiral*len)/(2*M_PI);
    zPeriod = stepSpiral;
    nPeriod = (_point.getZ() - zBegin)/zPeriod;
    if(fabs(zBegin + nPeriod*zPeriod - _point.getZ()) > fabs(zBegin + (nPeriod - 1)*zPeriod - _point.getZ())) nPeriod--;
    else if(fabs(zBegin + nPeriod*zPeriod - _point.getZ()) > fabs(zBegin + (nPeriod + 1)*zPeriod - _point.getZ())) nPeriod++;
    pointCenterCircle.setZ(zBegin + ((double)nPeriod)*zPeriod);
    /* Если центр окружности выходит за допустимые пределы, то точка не попадает в спираль */
    if(pointCenterCircle.getZ() > len/2.0 || pointCenterCircle.getZ() < -len/2.0) return false;
    /* Сравнение расстояния от центра окружности до определяемой точки и радиуса толщины спирали */
    rad = (_point.getX() - pointCenterCircle.getX())*(_point.getY() - pointCenterCircle.getY()) +
            (_point.getY() - pointCenterCircle.getY())*(_point.getY() - pointCenterCircle.getY()) +
            (_point.getZ() - pointCenterCircle.getZ())*(_point.getZ() - pointCenterCircle.getZ());
    if(rad > radDepth*radDepth*(1.0 + EPSILON_MUL)*(1.0 + EPSILON_MUL)) return false;
    else return true;
}

bool MathObject::isPointInFigure(const MathPoint3D &point)
{
    MathPoint3D translatePoint;
    double rad;

    translatePoint = translateBackPoint(point);
    switch(typeFigure)
    {
    case FIGURE_CILINDR:
        if(translatePoint.getZ() < -vectorParameters[1]/2.0 || translatePoint.getZ() > vectorParameters[1]/2.0) return false;
        rad = sqrt(translatePoint.getX()*translatePoint.getX() + translatePoint.getY()*translatePoint.getY());
        if(rad > vectorParameters[0]*(1.0 + EPSILON_MUL/1000.0)) return false;
        else return true;
        break;

    case FIGURE_SPHERE:
        rad = sqrt(translatePoint.getX()*translatePoint.getX() + translatePoint.getY()*translatePoint.getY() + translatePoint.getZ()*translatePoint.getZ());
        if(rad > vectorParameters[0]*(1.0 + EPSILON_MUL/1000.0)) return false;
        else return true;
        break;

    case FIGURE_SPIRAL:
        if(translatePoint.getZ() < -(vectorParameters[4]/2.0 + vectorParameters[0]) ||
                translatePoint.getZ() > (vectorParameters[4]/2.0 + vectorParameters[0])) return false;
        rad = sqrt(translatePoint.getX()*translatePoint.getX() + translatePoint.getY()*translatePoint.getY());
        if(rad > (vectorParameters[1] + vectorParameters[0])*(1.0 + EPSILON_MUL/1000.0)) return false;
        if(rad < (vectorParameters[1] - vectorParameters[0])*(1.0 - EPSILON_MUL/1000.0)) return false;
        if(isPointInSpiral(translatePoint)) return true;
        else return false;
        break;
    }
    return false;
}

MathPoint3D MathObject::translateBackPoint(const MathPoint3D &point)
{
    MathPoint3D translatePoint;
    uint32_t i, index;

    translatePoint = (MathPoint3D)point - centerFigure;
    for(i=0; i<vectorTypeRotates.size(); i++)
    {
        index = vectorTypeRotates.size() - i - 1;
        switch(vectorTypeRotates[index])
        {
        case AXIS_X:
            translatePoint = translatePoint.rotateParallelAxisX(MathPoint3D(0.0, 0.0, 0.0), -vectorValueRotates[index]);
            break;

        case AXIS_Y:
            translatePoint = translatePoint.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), -vectorValueRotates[index]);
            break;

        case AXIS_Z:
            translatePoint = translatePoint.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), -vectorValueRotates[index]);
            break;
        }
    }
    translatePoint = translatePoint.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), -sinNormal[1], cosNormal[1]);
    translatePoint = translatePoint.rotateParallelAxisY(MathPoint3D(0.0, 0.0, 0.0), -sinNormal[0], cosNormal[0]);
    return translatePoint;
}


//EOF
