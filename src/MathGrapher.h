/************************************************************************************************
 *    Класс выполняет отрисовку мат. объектов в трехмерном пространстве при помощи OpenGL       *
 ************************************************************************************************/
#ifndef MATHGRAPHER_H
#define MATHGRAPHER_H

#include <QtOpenGL/QGLWidget>
#include <QColor>
#include <vector>
#include "MathPoint3D.h"
#include "MathVector3D.h"
#include "MathGeomObjects.h"

#define ENABLE_LIGHTING
#define AXIS_LEN (0.2)


/* типы элементов для отрисовки и настройки - подставляются в MathGrapher::setEnable и MathGrapher::setDisable */
#define POINTS_TYPE     0
#define LINES_TYPE      1
#define VECTOR_TYPES    2
#define TETRAEDRS_TYPE  3
#define COORDINATE_GRID 4
#define AUTO_GRID       5
#define DRAW            6

/* Класс для отрисовки элементов в трехмерном пространстве */
class MathGrapher : public QGLWidget
{
    Q_OBJECT
public:
    /* Если использовать этот конструктор, то область и начало отсчета вычисляются автоматически */
    MathGrapher();

    /* Конструктор принимает значения начала отсчета и размеры области отрисовки */
    MathGrapher(const MathPoint3D& beginPoint, const double& dx, const double& dy, const double& dz);

    /* Настройка области отрисовки и начала отсчета */
    void setSizeRegion(const double& dx, const double& dy, const double& dz);
    void setBeginPoint(const MathPoint3D& point);

    /* Добавление точки(точек) и настройка их отображения */
    void addPoint(const MathPoint3D &point, const QColor& color);
    void addPoints(const MathPoint3D *points, uint32_t num, const QColor& color);
    void setColorPoints(const QColor& color);
    void clearPoints();

    /* Добавление линии(линий) и настройка их отображения */
    void addLine(const MathPoint3D &pointBegin, const MathPoint3D &pointEnd, const QColor &color);
    void addLines(const MathPoint3D *pointsBegin, const MathPoint3D *pointsEnd, const QColor &color, uint32_t num);
    void enableVolumeLines(bool yes);
    void setRadiusVolumeLines(double rad);
    void clearLines();

    /* Добавление вектора(векторов) и настройка их отображения */
    void addVector(const MathVector3D &vector, const QColor& color);                                /* Вектор с нормированным размером */
    void addVector(const MathVector3D &vector, const double& mulLenVector, const QColor& color);    /* Вектор с произвольным размером */
    void addVectors(const MathVector3D *vectors, uint32_t num, const QColor& color);
    double lenEightVector();                                                                        /* Длина единичного вектора */
    void clearVectors();

    /* Добавление тетраэдра */
    void addTetraedrs(const MathPoint3D& point1, const MathPoint3D& point2, const MathPoint3D& point3, const MathPoint3D& point4, const QColor& color);
    void clearTetraedrs();

    /* Включить/выключить */
    bool setEnable(uint32_t typeObject);
    bool setDisable(uint32_t typeObject);

    bool isEmpty();

    bool saveImage(const QString& fileName);
private:
    void setDefaultParameters();
    void calculateParameters(const MathPoint3D &point);
    /* параметры для OpenGL */
    double horAngle;
    double verAngle;
    double scale;
    double dxMoveRegion, dyMoveRegion, dzMoveRegion;
    int mouseCurX;
    int mouseCurY;
    double mulX;
    double mulY;

    /* общие параметры отрисовки */
    bool drawEnable;
    bool autoGridEnDraw;
    bool gridEnDraw;
    bool pointsEnDraw;
    bool linesEnDraw;
    bool vectorsEnDraw;
    bool tetraedrsEnDraw;
    bool enVolumeLines;

    MathPoint3D beginPointDraw;
    double dxDraw, dyDraw, dzDraw;
    bool allOnePoint;

    double radLines;

    /* массивы элементов */
    std::vector<double> vertexPointsDraw;
    std::vector<double> colorsPointsDraw;
    std::vector<double> normalsPointsDraw;
    std::vector<double> vertexLinesDraw;
    std::vector<double> colorsLinesDraw;
    std::vector<double> normalsLinesDraw;
    std::vector<double> vertexVectorsDraw;
    std::vector<double> colorsVectorsDraw;
    std::vector<double> normalsVectorsDraw;
    std::vector<double> vertexTetraedrs;
    std::vector<double> colorsTetraedrs;
    std::vector<double> normalsTetraedrs;
protected:
    virtual void postDraw();
    virtual void initializeGL();
    virtual void resizeGL(int width, int height);
    virtual void paintGL();
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent *event);
    virtual void keyPressEvent(QKeyEvent *event);

};

#endif // MATHGRAPHER_H
