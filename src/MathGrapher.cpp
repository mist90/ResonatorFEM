#include <stdint.h>
#include <GL/gl.h>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>
#include "MathGrapher.h"
#include "GLFunctions.h"

MathGrapher::MathGrapher()
{
    setDefaultParameters();

    setMinimumSize(640,480);
    setFocusPolicy(Qt::StrongFocus);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

MathGrapher::MathGrapher(const MathPoint3D& beginPoint, const double &dx, const double &dy, const double &dz)
{
    setDefaultParameters();
    setSizeRegion(dx, dy, dz);
    setBeginPoint(beginPoint);

    setMinimumSize(640,480);
    setFocusPolicy(Qt::StrongFocus);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void MathGrapher::setSizeRegion(const double &dx, const double &dy, const double &dz)
{
    autoGridEnDraw = false;
    dxDraw = dx;
    dyDraw = dy;
    dzDraw = dz;
    scale = 10.0/(dxDraw + dyDraw + dzDraw);
}

void MathGrapher::setBeginPoint(const MathPoint3D &point)
{
    autoGridEnDraw = false;
    beginPointDraw = point;
}

/* Добавление точки(точек) и настройка их отображения */
void MathGrapher::addPoint(const MathPoint3D &point, const QColor& color)
{
    uint32_t i, oldSize;
    setDisable(DRAW);
	calculateParameters(point);
    oldSize = vertexPointsDraw.size()/3;
    GetSphereVertex(radLines*3.0, point, vertexPointsDraw, normalsPointsDraw);
    for(i=0; i<(vertexPointsDraw.size()/3 - oldSize); i++)
    {
        colorsPointsDraw.push_back(color.red()/255.0);
        colorsPointsDraw.push_back(color.green()/255.0);
        colorsPointsDraw.push_back(color.blue()/255.0);
        colorsPointsDraw.push_back(((double)color.alpha())/255.0);
    }
}

void MathGrapher::addPoints(const MathPoint3D *points, uint32_t num, const QColor& color)
{
    uint32_t i;
    setDisable(DRAW);
    for(i=0; i<num; i++)
    {
	    calculateParameters(points[i]);
        addPoint(points[i], color);
    }
}

void MathGrapher::clearPoints()
{
    vertexPointsDraw.clear();
    colorsPointsDraw.clear();
    normalsPointsDraw.clear();

}

/* Добавление линии(линий) и настройка их отображения */
void MathGrapher::addLine(const MathPoint3D &pointBegin, const MathPoint3D &pointEnd, const QColor &color)
{
    uint32_t i, j;
    MathPoint3D point1, point2;
    MathVector3D vectorBegin, vectorEnd;
    double rad = radLines;

    setDisable(DRAW);
	calculateParameters(pointBegin);
    if(!enVolumeLines)
    {
        vertexLinesDraw.push_back(((MathPoint3D)pointBegin).getX());
        vertexLinesDraw.push_back(((MathPoint3D)pointBegin).getY());
        vertexLinesDraw.push_back(((MathPoint3D)pointBegin).getZ());
        vertexLinesDraw.push_back(((MathPoint3D)pointEnd).getX());
        vertexLinesDraw.push_back(((MathPoint3D)pointEnd).getY());
        vertexLinesDraw.push_back(((MathPoint3D)pointEnd).getZ());
        for(i=0; i<2; i++)
        {
            colorsLinesDraw.push_back(color.red()/255.0);
            colorsLinesDraw.push_back(color.green()/255.0);
            colorsLinesDraw.push_back(color.blue()/255.0);
            colorsLinesDraw.push_back(((double)color.alpha())/255.0);
        }
    }
    else
    {
        vectorBegin = MathVector3D((MathPoint3D)pointBegin, ((MathPoint3D)pointEnd));
        vectorEnd = vectorBegin;
        vectorEnd.setBegin(pointEnd);

        for(i=0; i<10; i++)
        {
            point1 = GetPointFromCircle(rad, vectorBegin, 36.0*i);
            vertexLinesDraw.push_back(point1.getX());
            vertexLinesDraw.push_back(point1.getY());
            vertexLinesDraw.push_back(point1.getZ());
            point2 = GetPointFromCircle(rad, vectorBegin, 36.0*(i + 1));
            vertexLinesDraw.push_back(point2.getX());
            vertexLinesDraw.push_back(point2.getY());
            vertexLinesDraw.push_back(point2.getZ());
            point2 = (point1 + point2)/2.0;
            for(j=0; j<4; j++)
            {
                normalsLinesDraw.push_back((point2.getX() - ((MathPoint3D)pointBegin).getX()));
                normalsLinesDraw.push_back((point2.getY() - ((MathPoint3D)pointBegin).getY()));
                normalsLinesDraw.push_back((point2.getZ() - ((MathPoint3D)pointBegin).getZ()));
            }
            point1 = GetPointFromCircle(rad, vectorEnd, 36.0*(i + 1));
            vertexLinesDraw.push_back(point1.getX());
            vertexLinesDraw.push_back(point1.getY());
            vertexLinesDraw.push_back(point1.getZ());
            point1 = GetPointFromCircle(rad, vectorEnd, 36.0*i);
            vertexLinesDraw.push_back(point1.getX());
            vertexLinesDraw.push_back(point1.getY());
            vertexLinesDraw.push_back(point1.getZ());
        }
        for(i=0; i<4*10; i++)
        {
            colorsLinesDraw.push_back(color.red()/255.0);
            colorsLinesDraw.push_back(color.green()/255.0);
            colorsLinesDraw.push_back(color.blue()/255.0);
            colorsLinesDraw.push_back(((double)color.alpha())/255.0);
        }
    }
    calculateParameters(pointEnd);
}

void MathGrapher::addLines(const MathPoint3D *pointsBegin, const MathPoint3D *pointsEnd, const QColor &color, uint32_t num)
{
    uint32_t i;
    for(i=0; i<num; i++) addLine(pointsBegin[i], pointsEnd[i], color);
}

void MathGrapher::enableVolumeLines(bool yes)
{
    enVolumeLines = yes;
    clearLines();
}

void MathGrapher::setRadiusVolumeLines(double rad)
{
    radLines = rad;
}

void MathGrapher::clearLines()
{
    vertexLinesDraw.clear();
    colorsLinesDraw.clear();
    normalsLinesDraw.clear();
}

/* Добавление вектора(векторов) и настройка их отображения */
void MathGrapher::addVector(const MathVector3D &vector, const QColor& color)
{
    addVector((MathVector3D(vector)).unitVector(), radLines*20.0, color);
}

void MathGrapher::addVector(const MathVector3D &vector, const double &mulLenVector, const QColor &color)
{
    MathPoint3D begin, end;
    MathVector3D vectorFixLen;
    uint32_t i;

    setDisable(DRAW);
    vectorFixLen = ((MathVector3D)vector)*mulLenVector;
    calculateParameters(vectorFixLen.getBegin());
    GetConusVertex(mulLenVector/20.0, vectorFixLen.getBegin(), vectorFixLen.getEnd(), vertexVectorsDraw, normalsVectorsDraw, 15);
    for(i=0; i<15*3; i++)
    {
        colorsVectorsDraw.push_back(color.red()/255.0);
        colorsVectorsDraw.push_back(color.green()/255.0);
        colorsVectorsDraw.push_back(color.blue()/255.0);
        colorsVectorsDraw.push_back(((double)color.alpha())/255.0);
    }
    calculateParameters(vectorFixLen.getEnd());
    end = vectorFixLen.getEnd();
    vectorFixLen = vectorFixLen*(0.6);
    begin = vectorFixLen.getEnd();
    GetConusVertex(mulLenVector*3.0/20.0, begin, end, vertexVectorsDraw, normalsVectorsDraw, 15);
    for(i=0; i<15*3; i++)
    {
        colorsVectorsDraw.push_back(color.red()/255.0);
        colorsVectorsDraw.push_back(color.green()/255.0);
        colorsVectorsDraw.push_back(color.blue()/255.0);
        colorsVectorsDraw.push_back(((double)color.alpha())/255.0);
    }
}

void MathGrapher::addVectors(const MathVector3D *vectors, uint32_t num, const QColor& color)
{
    uint32_t i;
    MathPoint3D begin, end;
    setDisable(DRAW);
    for(i=0; i<num; i++)
    {
	    begin = ((MathVector3D)vectors[i]).getBegin();
        end = ((MathVector3D)vectors[i]).getEnd();
        calculateParameters(begin);
        addVector(vectors[i], color);
        calculateParameters(end);
    }
}

double MathGrapher::lenEightVector()
{
    return radLines*30.0;
}

void MathGrapher::addTetraedrs(const MathPoint3D &point1, const MathPoint3D &point2, const MathPoint3D &point3, const MathPoint3D &point4, const QColor &color)
{
    MathPoint3D *points[4];
    MathVector3D vec1, vec2, vecNormal;
    uint32_t i, j;

    setDisable(DRAW);
    points[0] = (MathPoint3D*)&point1;
    points[1] = (MathPoint3D*)&point2;
    points[2] = (MathPoint3D*)&point3;
    points[3] = (MathPoint3D*)&point4;

    for(i=0; i<4; i++)
    {
        for(j=0; j<3; j++)
        {
            vertexTetraedrs.push_back(points[(i + j)%4]->getX());
            vertexTetraedrs.push_back(points[(i + j)%4]->getY());
            vertexTetraedrs.push_back(points[(i + j)%4]->getZ());
        }
        vec1 = MathVector3D(*points[i%4], *points[(i + 1)%4]);
        vec2 = MathVector3D(*points[i%4], *points[(i + 2)%4]);
        vecNormal = vec1^vec2;
        vecNormal = vecNormal.unitVector();
        vec1 = MathVector3D(*points[i%4], *points[(i + 3)%4]);
        if(vec1*vecNormal > 0.0) vecNormal = MathVector3D(0.0, 0.0, 0.0) - vecNormal;
        for(j=0; j<3; j++)
        {
            normalsTetraedrs.push_back(vecNormal.getX());
            normalsTetraedrs.push_back(vecNormal.getY());
            normalsTetraedrs.push_back(vecNormal.getZ());
        }
        calculateParameters(*points[i]);
    }

    for(i=0; i<12; i++)
    {
        colorsTetraedrs.push_back(color.red());
        colorsTetraedrs.push_back(color.green());
        colorsTetraedrs.push_back(color.blue());
        colorsTetraedrs.push_back(((double)color.alpha())/255.0);
    }

}

void MathGrapher::clearTetraedrs()
{
    vertexTetraedrs.clear();
    colorsTetraedrs.clear();
    normalsTetraedrs.clear();
}

void MathGrapher::clearVectors()
{
    vertexVectorsDraw.clear();
    colorsVectorsDraw.clear();
    normalsVectorsDraw.clear();
}

bool MathGrapher::setEnable(uint32_t typeObject)
{
    switch(typeObject)
    {
    case POINTS_TYPE:
        pointsEnDraw = true;
        break;

    case LINES_TYPE:
        linesEnDraw = true;
        break;

    case VECTOR_TYPES:
        vectorsEnDraw = true;
        break;

    case TETRAEDRS_TYPE:
        tetraedrsEnDraw = true;

    case COORDINATE_GRID:
        gridEnDraw = true;
        break;

    case AUTO_GRID:
        autoGridEnDraw = true;
        break;

    case DRAW:
        drawEnable = true;
        updateGL();
        break;

    default:
        return false;
    }
    return true;
}

bool MathGrapher::setDisable(uint32_t typeObject)
{
    switch(typeObject)
    {
    case POINTS_TYPE:
        pointsEnDraw = false;
        break;

    case LINES_TYPE:
        linesEnDraw = false;
        break;

    case VECTOR_TYPES:
        vectorsEnDraw = false;
        break;

    case TETRAEDRS_TYPE:
        tetraedrsEnDraw = false;

    case COORDINATE_GRID:
        gridEnDraw = false;
        break;

    case AUTO_GRID:
        autoGridEnDraw = false;
        break;

    case DRAW:
        drawEnable = false;
        break;

    default:
        return false;
    }
    return true;
}

bool MathGrapher::isEmpty()
{
    if(!vertexPointsDraw.size() && !vertexLinesDraw.size() && !vertexVectorsDraw.size() && !vertexTetraedrs.size()) return true;
    else return false;
}

bool MathGrapher::saveImage(const QString &fileName)
{
    QImage image = grabFrameBuffer();

    return image.save(fileName);
}

void MathGrapher::setDefaultParameters()
{
    horAngle = 0.0;
    verAngle = -90.0;
    scale = 1.0;
    dxMoveRegion = dyMoveRegion = 0.0;
    dzMoveRegion = -50.0;
    mulX = mulY = 1.0;

    autoGridEnDraw = gridEnDraw = pointsEnDraw = linesEnDraw = vectorsEnDraw = true;
    enVolumeLines = true;
    drawEnable = false;
    radLines = 0.01;
    allOnePoint = false;
    setSizeRegion(1.0, 1.0, 1.0);
    setBeginPoint(MathPoint3D(0.0, 0.0, 0.0));
    setEnable(AUTO_GRID);
    setEnable(POINTS_TYPE);
    setEnable(LINES_TYPE);
    setEnable(VECTOR_TYPES);
    setEnable(TETRAEDRS_TYPE);
    setEnable(COORDINATE_GRID);
    setEnable(AUTO_GRID);
}

void MathGrapher::calculateParameters(const MathPoint3D &point)
{
    MathPoint3D _point = point;
    if(!autoGridEnDraw) return;
    if(isEmpty())
    {
        dxDraw = dyDraw = dzDraw = 1.0;
        beginPointDraw.setX(_point.getX() - dxDraw);
        beginPointDraw.setY(_point.getY() - dyDraw);
        beginPointDraw.setZ(_point.getZ() - dzDraw);
        scale = 10.0/(dxDraw + dyDraw + dzDraw);
        allOnePoint = true;
        return;
    }
    if(allOnePoint)
    {
        beginPointDraw.setX(beginPointDraw.getX() + dxDraw);
        beginPointDraw.setY(beginPointDraw.getY() + dyDraw);
        beginPointDraw.setZ(beginPointDraw.getZ() + dzDraw);
        dxDraw = dyDraw = dzDraw = 0.0;
        allOnePoint = false;
    }
    if(_point.getX() >= beginPointDraw.getX() && _point.getX() <= beginPointDraw.getX() + dxDraw &&
       _point.getY() >= beginPointDraw.getY() && _point.getY() <= beginPointDraw.getY() + dyDraw &&
       _point.getZ() >= beginPointDraw.getZ() && _point.getZ() <= beginPointDraw.getZ() + dzDraw)
        return;
    /* если новая точка не лежит в выделенной области */
    /* по X */
    if(_point.getX() < beginPointDraw.getX())
    {
        beginPointDraw.setX(_point.getX());
        dxDraw += beginPointDraw.getX() - _point.getX();
    }
    else if(_point.getX() > beginPointDraw.getX() + dxDraw)
    {
        dxDraw = _point.getX() - beginPointDraw.getX();
    }
    /* по Y */
    if(_point.getY() < beginPointDraw.getY())
    {
        beginPointDraw.setY(_point.getY());
        dyDraw += beginPointDraw.getY() - _point.getY();
    }
    else if(_point.getY() > beginPointDraw.getY() + dyDraw)
    {
        dyDraw = _point.getY() - beginPointDraw.getY();
    }
    /* по Z */
    if(_point.getZ() < beginPointDraw.getZ())
    {
        beginPointDraw.setZ(_point.getZ());
        dzDraw += beginPointDraw.getZ() - _point.getZ();
    }
    else if(_point.getZ() > beginPointDraw.getZ() + dzDraw)
    {
        dzDraw = _point.getZ() - beginPointDraw.getZ();
    }

    scale = 10.0/(dxDraw + dyDraw + dzDraw);
}

void MathGrapher::postDraw()
{

}

void MathGrapher::initializeGL()
{
    GLfloat lightAmbient[]={0.0f,0.0f,0.0f,1.0f};
    GLfloat lightDiffuse[]={0.5f,0.5f,0.5f,1.0f};
    GLfloat lightSpec[]={0.0f,0.0f,0.0f,1.0f};

    GLfloat materialAmbient[] = {0.0f, 0.0f, 0.0f, 1.0f};
    GLfloat materialSpecular[] = {0.0f, 0.0f, 0.0f, 1.0f};
    GLfloat shinines = 250.0;

    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_ALPHA_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_NORMALIZE);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,materialAmbient);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,materialSpecular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shinines);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

    glLightfv(GL_LIGHT0,GL_AMBIENT,lightAmbient);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,lightDiffuse);
    glLightfv(GL_LIGHT0,GL_SPECULAR,lightSpec);

    glEnable(GL_LIGHT0);
}

void MathGrapher::resizeGL(int width, int height)
{
    if(width > height)
        mulX = (double)(width) / (double)(height);
    else mulY = (double)(height) / (double)(width);
    glViewport(0, 0, (GLint)width, (GLint)height);
    updateGL();
}

void MathGrapher::paintGL()
{
    uint32_t optSize;
    double sizeFont = AXIS_LEN/10.0;
    double dxCurrent, dyCurrent, dzCurrent;
    GLfloat lightPos[4] = {1.0, 1.0, 1.0, 0.0};

    optSize = dxDraw;
    if(optSize > dyDraw) optSize = dyDraw;
    if(optSize > dzDraw) optSize = dzDraw;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1.0*mulX, 1.0*mulX, -1.0*mulY, 1.0*mulY, 10.0, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if(drawEnable)
    {
        dxCurrent = -(beginPointDraw.getX() + dxDraw/2.0);
        dyCurrent = -(beginPointDraw.getY() + dyDraw/2.0);
        dzCurrent = -(beginPointDraw.getZ() + dzDraw/2.0);
        glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
        glTranslated(dxMoveRegion, dyMoveRegion, dzMoveRegion);
        glScaled(scale, scale, scale);
        glRotated(verAngle, 1.0, 0.0, 0.0);
        glRotated(horAngle, 0.0, 0.0, 1.0);
        glTranslated(dxCurrent, dyCurrent, dzCurrent);
        /* рисование координатных осей */
        glPushMatrix();
        glLoadIdentity();
        glTranslated(0, 0.0, -12.0);
        glTranslated(-mulX, -mulY, 0.0);
        glRotated(verAngle, 1.0, 0.0, 0.0);
        glRotated(horAngle, 0.0, 0.0, 1.0);

        glBegin(GL_LINES);
        glColor3d(0.0, 0.0, 0.0);
        /* Z */
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(0.0, 0.0, AXIS_LEN);
        /* Y */
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(0.0, AXIS_LEN, 0.0);
        /* X */
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(AXIS_LEN, 0.0, 0.0);

        /* рисование подписй к осям */
        glColor3d(0.0, 0.0, 255.0);
        /* Z */
        glVertex3d(0.0 - sizeFont, 0.0, AXIS_LEN + 4*sizeFont);
        glVertex3d(0.0 + sizeFont, 0.0, AXIS_LEN + 4*sizeFont);
        glVertex3d(0.0 + sizeFont, 0.0, AXIS_LEN + 4*sizeFont);
        glVertex3d(0.0 - sizeFont, 0.0, AXIS_LEN);
        glVertex3d(0.0 - sizeFont, 0.0, AXIS_LEN);
        glVertex3d(0.0 + sizeFont, 0.0, AXIS_LEN);
        /* Y */
        glVertex3d(0.0, AXIS_LEN - sizeFont, 0.0);
        glVertex3d(0.0, AXIS_LEN, - 2*sizeFont);
        glVertex3d(0.0, AXIS_LEN + sizeFont, 0.0);
        glVertex3d(0.0, AXIS_LEN - sizeFont, 0.0 - 4*sizeFont);
        /* X */
        glVertex3d(AXIS_LEN - sizeFont, 0.0, 0.0);
        glVertex3d(AXIS_LEN + sizeFont, 0.0, 0.0 - 4*sizeFont);
        glVertex3d(AXIS_LEN + sizeFont, 0.0, 0.0);
        glVertex3d(AXIS_LEN - sizeFont, 0.0, 0.0 - 4*sizeFont);
        glEnd();
        glPopMatrix();
        /* рисование точек */
        if(pointsEnDraw)
        {
            glEnable(GL_LIGHTING);
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glEnableClientState((GL_NORMAL_ARRAY));
            glVertexPointer(3, GL_DOUBLE, 0, vertexPointsDraw.data());
            glColorPointer(4, GL_DOUBLE, 0, colorsPointsDraw.data());
            glNormalPointer(GL_DOUBLE, 0, normalsPointsDraw.data());
            glDrawArrays(GL_QUADS, 0, vertexPointsDraw.size()/3);
            glDisableClientState(GL_NORMAL_ARRAY);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisable(GL_LIGHTING);
        }

        /* рисование линий */
        if(linesEnDraw)
        {
            if(enVolumeLines) glEnable(GL_LIGHTING);
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glEnableClientState((GL_NORMAL_ARRAY));
            glVertexPointer(3, GL_DOUBLE, 0, vertexLinesDraw.data());
            glColorPointer(4, GL_DOUBLE, 0, colorsLinesDraw.data());
            glNormalPointer(GL_DOUBLE, 0, normalsLinesDraw.data());
            if(enVolumeLines) glDrawArrays(GL_QUADS, 0, vertexLinesDraw.size()/3);
            else glDrawArrays(GL_LINES, 0, vertexLinesDraw.size()/3);
            glDisableClientState(GL_NORMAL_ARRAY);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            if(enVolumeLines) glDisable(GL_LIGHTING);
        }

        /* рисование векторов */
        if(vectorsEnDraw)
        {
            glEnable(GL_LIGHTING);
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glVertexPointer(3, GL_DOUBLE, 0, vertexVectorsDraw.data());
            glColorPointer(4, GL_DOUBLE, 0, colorsVectorsDraw.data());
            glNormalPointer(GL_DOUBLE, 0, normalsVectorsDraw.data());
            glDrawArrays(GL_TRIANGLES, 0, vertexVectorsDraw.size()/3);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            glDisable(GL_LIGHTING);
        }

        /* рисование тетраэдров */
        if(tetraedrsEnDraw)
        {
            glEnable(GL_LIGHTING);
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glVertexPointer(3, GL_DOUBLE, 0, vertexTetraedrs.data());
            glColorPointer(4, GL_DOUBLE, 0, colorsTetraedrs.data());
            glNormalPointer(GL_DOUBLE, 0, normalsTetraedrs.data());
            glDrawArrays(GL_TRIANGLES, 0, vertexTetraedrs.size()/3);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            glDisable(GL_LIGHTING);
        }

    }
}

void MathGrapher::mousePressEvent(QMouseEvent *event)
{
    mouseCurX = event->pos().x();
    mouseCurY = event->pos().y();
}

void MathGrapher::mouseMoveEvent(QMouseEvent *event)
{
    int dx, dy;
    dx = event->pos().x() - mouseCurX;
    dy = event->pos().y() - mouseCurY;
    mouseCurX = event->pos().x();
    mouseCurY = event->pos().y();
    horAngle = horAngle + (double)dx * 3.14/18.0;
    verAngle = verAngle + (double)dy * 3.14/18.0;
    updateGL();
}

void MathGrapher::wheelEvent(QWheelEvent *event)
{
    scale *= pow(10, (double)event->delta()/800.0);
    updateGL();
}

void MathGrapher::keyPressEvent(QKeyEvent *event)
{
    switch(event->key())
    {
        case Qt::Key_Up:
        dyMoveRegion += 1.0;
        break;

        case Qt::Key_Down:
        dyMoveRegion -= 1.0;
        break;

        case Qt::Key_Left:
        dxMoveRegion -= 1.0;
        break;

        case Qt::Key_Right:
        dxMoveRegion += 1.0;
        break;
    }


    updateGL();
}

//EOF
