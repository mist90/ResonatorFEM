#ifndef MICROTETRAEDR_H
#define MICROTETRAEDR_H

#include "MathPoint3D.h"
#include "MathVector3D.h"
#include "MathMatrix.h"
#include "MicroNode.h"
#include "MicroFlags.h"
#include <stdint.h>
#include <list>
#include <vector>

#define EPSILON_RAD_QUAD (0.0000000001) /* погрешность вычисления сферы */
#define CALC_SPHERE_METHOD_SLAU         /* метод вычисления уравнения описанной сферы */

/* Класс, хранящий информацию о ребре тетраэдра */
class MicroEdge
{
public:
    MicroEdge();
    MicroEdge(const std::list<MicroNode>::iterator& node1, const std::list<MicroNode>::iterator& node2);
    void                            setNode1(const std::list<MicroNode>::iterator& node1);
    void                            setNode2(const std::list<MicroNode>::iterator& node2);
    std::list<MicroNode>::iterator  getNode1();
    std::list<MicroNode>::iterator  getNode2();
    /* Возвращает true, если два ребра опираются на одинаковый узлы */
    bool                            operator==(const MicroEdge& edge);
    /* Возвращает true, если два ребра опираются на одинаковые узлы крест-накрест;
       предварительно нужно проверить равны ли ребра вообще */
    bool                            isReverse(const MicroEdge& edge);
    /* Замена номеров узлов */
    void                            reverseEdge();
    double                          lenEdge();
    bool                            isFlags(uint32_t flags);
private:
    std::list<MicroNode>::iterator _node1;
    std::list<MicroNode>::iterator _node2;
};

class MicroTetraedr;
/* Контейнер для тетраэдров */
class ListMicroTetraedr:public std::list<MicroTetraedr>
{
public:
    ListMicroTetraedr();
    /* поиск тетраэдра, в который попала точка */
    bool findTetraedr(MathPoint3D &point, ListMicroTetraedr::iterator &it);
};

/* Структура для хранения тетраэдра и информации о смежных тетраэдрах */
class MicroTetraedr
{
public:
    MicroTetraedr();
    std::list<MicroNode>::iterator& nodes(uint32_t i);
    ListMicroTetraedr::iterator&    tetraedrs(uint32_t i);
    void                            clearTetraedrs();
    /* Принадлежит ли точка тетраэдру? */
    bool                            isInTetraedr(const MathPoint3D &point);
    /* Принадлежит ли точка той половине пространства, где находится тетраэдр, отн. грани numSurface */
    bool                            isInTetraedr(const MathPoint3D &point, uint32_t numSurface);
    /* Принадлежит ли точка узлу плоскости (треугольнику)? */
    bool                            isPointInSurface(const MathPoint3D &point, uint32_t numSurface);
    /* Принадлежит ли ребро границе плоскости (треугольнику)? */
    bool                            isEdgeInSurface(const MathPoint3D &point1, const MathPoint3D &point2, uint32_t numSurface);
    /* Являются ли точка узлом тетраэдра? */
    bool                            isNodeTetraedr(const MathPoint3D &point);
    /* Являются ли 3 точки узлами тетраэдра? */
    bool                            isNode3Tetraedr(const MathPoint3D &point1, const MathPoint3D &point2, const MathPoint3D &point3);
    /* Получение номера плоскости по трем точкам */
    bool                            getNumSurface(const MathPoint3D &point1, const MathPoint3D &point2, const MathPoint3D &point3, uint32_t& getNum);
    /* Получение номера плоскости по итератору на соседний тетраэдр */
    bool                            getNumSurface(ListMicroTetraedr::iterator it, uint32_t& getNum);
    /* Получение количества свободных граней */
    uint32_t                        getCountEmptyNeighbours();
    /* Сфера, описанная вокруг тетраэдра */
    void                            makeCircumSphere();
    /* Возвращает центр описанной сферы */
    MathPoint3D                     pointCenterSphere();
    /* Определение граничных граней - вызывать только после полного заполнения тетраэдра */
    void                            calculateBoundariesSurface();
    bool                            isBoundarySurface(uint32_t index);
    /* Флаги */
    bool                            isFlags(uint32_t numSurface, uint32_t flags);
    void                            addFlags(uint32_t numSurface, uint32_t flags);
    void                            clearFlags(uint32_t numSurface, uint32_t flags);
    /* Точка входит в сферу? */
    bool                            isInSphere(const MathPoint3D &point);
    /* Объем тетраэдра */
    double                          volume();
    /* Задание номера тетраэдра */
    void                            setSerialNumber(uint32_t number);
    uint32_t                        getSerialNumber();
                            /* Функции для МКЭ */
    /* Получение ребра тетраэдра по номеру  */
    uint32_t                        getNumBeginEdge(uint32_t numEdge);
    uint32_t                        getNumEndEdge(uint32_t numEdge);
    MicroEdge                       getEdge(uint32_t numEdge);
    uint32_t                        getNumEdgeFromSurface(uint32_t numSurface, uint32_t numEdge);
    uint32_t                        getNumPointFromSurface(uint32_t numSurface, uint32_t numPoint);
    /* Вычислить значение базисной функции numFunc в точке point */
    MathVector3D                    getValueBasisFunc(const MathPoint3D& point, uint32_t numEdge);
    /* Получить вектор, нормальный к поверхности numSurface */
    MathVector3D                    getNormVector(uint32_t numSurface);

private:
    std::list<MicroNode>::iterator      _nodes[4];
    ListMicroTetraedr::iterator         _tetraedrs[4];
    uint32_t                            _flags[4];
    MathPoint3D                         _pointCenterSphere;
    double                              _radSphereQuad;
    uint32_t                            _serialNumber;
};

class MicroTriangle
{
public:
    MicroTriangle(std::list<MicroNode>::iterator node1, std::list<MicroNode>::iterator node2, std::list<MicroNode>::iterator node3);
    std::list<MicroNode>::iterator   node(uint32_t index);
    bool                             isNodeTriangle(MathPoint3D &point);
    bool                             operator==(const MicroTriangle& triangle);
private:
    std::list<MicroNode>::iterator   _nodes[3];
};

/* Получение пустого итератора */
ListMicroTetraedr::iterator EmptyIterator();

/* Два тетраэдра совмещены? */
bool IsSharedTetraedrs(MicroTetraedr& tetraedr1, MicroTetraedr& tetraedr2);

/* Два тетраэдра пересекаются? */
bool IsCrossingTetraedrs(MicroTetraedr& tetraedr1, MicroTetraedr& tetraedr2);

/* Тетраэдры удовлетворяют условию Делоне? */
bool IsDeloneTetraedrs(ListMicroTetraedr::iterator &tetraedr1, ListMicroTetraedr::iterator &tetraedr2);

/* Тетраэдр удовлетворяют условию Делоне по отношению к соседям? */
bool IsDeloneTetraedr(ListMicroTetraedr::iterator &tetraedr);

/* Тетраэдр вырожден, epsilonVolume - минимальный объем */
bool IsSingularTetraedr(MathPoint3D &point1, MathPoint3D &point2, MathPoint3D &point3, MathPoint3D &point4, double epsilonVolume);

double VolumeTetraedr(MathPoint3D &point1, MathPoint3D &point2, MathPoint3D &point3, MathPoint3D &point4);

#endif // MICROTETRAEDR_H
