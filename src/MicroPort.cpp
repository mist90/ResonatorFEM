#include "MicroPort.h"

MicroPort::MicroPort()
{
    clear();
}

void MicroPort::makeCoaxial(const MathVector3D& base, const double& lenCilindr, const double& radMin, const double& radMax, const double step)
{
    MathVector3D vector(base);
    double deltaRad = radMax - radMin;
    double _step;

    clear();
    if(deltaRad < 0.0) return;
    if(step != AUTO_STEP && step < 0) return;
    _step = step;
    vector = (vector.unitVector())*lenCilindr;
    typePort = TYPE_PORT_COAXIAL;
    objectsPort.resize(5, MathObject());
    objectsPort[0].makeCilindr((vector*0.5).getEnd(), radMax, lenCilindr, _step);
    objectsPort[0].setVector(vector);
    objectsPort[1].makeCilindr((vector*0.5).getEnd(), radMax, lenCilindr, _step);
    objectsPort[1].setVector(vector);
    fieldSurface.makeCilindr(vector.getEnd(), radMax, lenCilindr/20.0, _step);
    fieldSurface.setVector(vector);
    objectsPort[2].makeCilindr((vector*0.5).getEnd(), radMin + deltaRad*1.0/3.0, lenCilindr, _step);
    objectsPort[2].setVector(vector);
    objectsPort[3].makeCilindr((vector*0.5).getEnd(), radMin + deltaRad*2.0/3.0, lenCilindr, _step);
    objectsPort[3].setVector(vector);
    if(step == AUTO_STEP)
        _step = 2.0;
    else _step = step;
    objectsPort[4].makeCilindr((vector*0.5).getEnd(), radMin, lenCilindr, _step);
    objectsPort[4].setVector(vector);
}

void MicroPort::setI(const double& maxAmpl, const double& phase)
{
    maxAmplitude = maxAmpl;
    beginPhase = phase;
}

void MicroPort::setRadiation(const bool &enable)
{
    isRadiation = enable;
}

bool MicroPort::getObjects(std::vector<MathObject>& objects, std::vector<double>& sigma, std::vector<double>& epsilon)
{
    uint32_t i;

    switch (typePort)
    {
    case TYPE_PORT_COAXIAL:
        for(i=0; i<objectsPort.size(); i++)
        {
            objects.push_back(objectsPort[i]);
            if(i == 0 || i == 4) sigma.push_back(SIGMA_IDEAL_METALL);
            else sigma.push_back(0.0);
            epsilon.push_back(1.0);
        }

        break;
    default:
        return false;
        break;
    }

    return true;
}

bool MicroPort::isPointInPort(const MathPoint3D &point)
{
    if(!objectsPort.size()) return false;
    switch(typePort)
    {
    case TYPE_PORT_COAXIAL:
        if(objectsPort[0].isPointInFigure(point)) return true;
        else return false;
    }

    return false;
}

bool MicroPort::isPointInEndPort(const MathPoint3D &point)
{
    if(!objectsPort.size()) return false;
    if(fieldSurface.isPointInFigure(point)) return true;
    else return false;
}

bool MicroPort::getFieldPoints(const MathPoint3D &point, MathVector3D &amplVectorH, double &phaseVectorH)
{
    MathPoint3D translatePoint, endPoint;

    double rad;

    if(!objectsPort.size()) return false;
    if(!isRadiation) return false;
    if(fieldSurface.isPointInFigure(point))
    {
        switch (typePort) {
        case TYPE_PORT_COAXIAL:
            translatePoint = fieldSurface.translateBackPoint(point);
            translatePoint.setZ(0.0);
            rad = sqrt(translatePoint.getX()*translatePoint.getX() + translatePoint.getY()*translatePoint.getY());
            endPoint = (translatePoint.rotateParallelAxisZ(MathPoint3D(0.0, 0.0, 0.0), 90.0)) + translatePoint;
            amplVectorH = MathVector3D(translatePoint, endPoint);
            amplVectorH = amplVectorH.unitVector();
            amplVectorH = amplVectorH*(maxAmplitude/(2.0*M_PI*rad));
            amplVectorH.setBegin(point);
            phaseVectorH = beginPhase;
            break;
        default:
            return false;
        }
        return true;
    }
    else return false;
}

bool MicroPort::getAbsorbPoints(const MathPoint3D &point)
{
    if(!objectsPort.size()) return false;
    if(isRadiation) return false;
    if(fieldSurface.isPointInFigure(point)) return true;
    else return false;
}

void MicroPort::clear()
{
    maxAmplitude = 1.0;
    beginPhase = 0.0;
    typePort = TYPE_PORT_COAXIAL;
    isRadiation = true;
    objectsPort.clear();
    fieldSurface.clear();
}


//EOF
