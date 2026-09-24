#ifndef MICROPORT
#define MICROPORT

#include <vector>
#include <list>
#include <stdint.h>
#include "MathObject.h"
#include "MicroTetraedr.h"

/* Форма порта */
#define TYPE_PORT_COAXIAL   (0)
/* Идеальный металл */
#define SIGMA_IDEAL_METALL (-1)


class MicroPort
{
public:
    MicroPort();
    void makeCoaxial(const MathVector3D &base, const double &lenCilindr, const double &radMin, const double &radMax, const double step = AUTO_STEP);
    void setI(const double &maxAmpl, const double &phase);
    void setRadiation(const bool &enable);
    bool getObjects(std::vector<MathObject>& objects, std::vector<double>& sigma, std::vector<double>& epsilon);
    bool isPointInPort(const MathPoint3D& point);
    bool isPointInEndPort(const MathPoint3D& point);
    bool getFieldPoints(const MathPoint3D &point, MathVector3D &amplVectorH, double &phaseVectorH);
    bool getAbsorbPoints(const MathPoint3D& point);
    void clear();
private:
    uint32_t typePort;                          /* тип порта, например коаксиальный */
    bool isRadiation;                           /* излучающий/поглощающий */
    std::vector<MathObject> objectsPort;        /* геометрические объекты порта */
    MathObject fieldSurface;
    double maxAmplitude;                        /* максимальное значение амплитуды поля, В/м */
    double beginPhase;                          /* начальная фаза поля, рад */
};


#endif // MICROPORT

