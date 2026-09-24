#include "GLFunctions.h"
#include <math.h>

void GetSphereVertex(GLdouble r, MathPoint3D point, std::vector<double>& points, std::vector<double> &normals, unsigned int numStepVer, unsigned int numStepRad)
{
    double dAngle1 = 180.0/(double)numStepVer, dAngle2 = 360.0/(double)numStepRad, sumAngle1, sumAngle2, pi=(double)M_PI;
    unsigned int i, j;

    sumAngle2 = 0;
    for(i=0; i<numStepRad; i++)
    {
        sumAngle1 = -90;
        for(j=0; j<numStepVer; j++)
        {
            points.push_back(r*cos(sumAngle1*pi/180)*cos(sumAngle2*pi/180) + point.getX());
            points.push_back(r*sin(sumAngle1*pi/180) + point.getY());
            points.push_back(r*cos(sumAngle1*pi/180)*sin(sumAngle2*pi/180) + point.getZ());
            normals.push_back((points[points.size() - 3] - point.getX())/r);
            normals.push_back((points[points.size() - 2] - point.getY())/r);
            normals.push_back((points[points.size() - 1] - point.getZ())/r);

            points.push_back(r*cos((sumAngle1 + dAngle1)*pi/180)*cos(sumAngle2*pi/180) + point.getX());
            points.push_back(r*sin((sumAngle1 + dAngle1)*pi/180) + point.getY());
            points.push_back(r*cos((sumAngle1 + dAngle1)*pi/180)*sin(sumAngle2*pi/180) + point.getZ());
            normals.push_back((points[points.size() - 3] - point.getX())/r);
            normals.push_back((points[points.size() - 2] - point.getY())/r);
            normals.push_back((points[points.size() - 1] - point.getZ())/r);

            points.push_back(r*cos((sumAngle1 + dAngle1)*pi/180)*cos((sumAngle2+dAngle2)*pi/180) + point.getX());
            points.push_back(r*sin((sumAngle1 + dAngle1)*pi/180) + point.getY());
            points.push_back(r*cos((sumAngle1 + dAngle1)*pi/180)*sin((sumAngle2+dAngle2)*pi/180) + point.getZ());
            normals.push_back((points[points.size() - 3] - point.getX())/r);
            normals.push_back((points[points.size() - 2] - point.getY())/r);
            normals.push_back((points[points.size() - 1] - point.getZ())/r);

            points.push_back(r*cos(sumAngle1*pi/180)*cos((sumAngle2+dAngle2)*pi/180) + point.getX());
            points.push_back(r*sin(sumAngle1*pi/180) + point.getY());
            points.push_back(r*cos(sumAngle1*pi/180)*sin((sumAngle2+dAngle2)*pi/180) + point.getZ());
            normals.push_back((points[points.size() - 3] - point.getX())/r);
            normals.push_back((points[points.size() - 2] - point.getY())/r);
            normals.push_back((points[points.size() - 1] - point.getZ())/r);

            sumAngle1 = sumAngle1 + dAngle1;
        }
        sumAngle2 = sumAngle2 + dAngle2;
    }

}

void DrawSphereTexture(GLdouble r, MathPoint3D point, GLuint numTextures, GLdouble xozAngle, unsigned int numStepVer,unsigned int numStepRad)
{
    GLdouble *rCoord;
    GLdouble *tCoord;
    GLdouble *nCoord;
    glBindTexture(GL_TEXTURE_2D, numTextures);
    double dAngle1 = 180/(double)numStepVer, dAngle2 = 360/(double)numStepRad, sumAngle1, sumAngle2, pi=(double)M_PI;
    unsigned int i, j, nRect, in;
    nRect = numStepVer*numStepRad;
    rCoord = new GLdouble[nRect*12];
    tCoord = new GLdouble[nRect*12];
    nCoord = new GLdouble[nRect*12];
    sumAngle2 = xozAngle;
    in = 0;
    for(i=0; i<numStepRad; i++)
    {
        sumAngle1 = -90;
        for(j=0; j<numStepVer; j++)
        {
            rCoord[in] = r*cos(sumAngle1*pi/180)*cos(sumAngle2*pi/180) + point.getX();
            rCoord[in + 1] = r*sin(sumAngle1*pi/180) + point.getY();
            rCoord[in + 2] = r*cos(sumAngle1*pi/180)*sin(sumAngle2*pi/180) + point.getZ();

            rCoord[in + 3] = r*cos((sumAngle1 + dAngle1)*pi/180)*cos(sumAngle2*pi/180) + point.getX();
            rCoord[in + 4] = r*sin((sumAngle1 + dAngle1)*pi/180) + point.getY();
            rCoord[in + 5] = r*cos((sumAngle1 + dAngle1)*pi/180)*sin(sumAngle2*pi/180) + point.getZ();

            rCoord[in + 6] = r*cos((sumAngle1 + dAngle1)*pi/180)*cos((sumAngle2+dAngle2)*pi/180) + point.getX();
            rCoord[in + 7] = r*sin((sumAngle1 + dAngle1)*pi/180) + point.getY();
            rCoord[in + 8] = r*cos((sumAngle1 + dAngle1)*pi/180)*sin((sumAngle2+dAngle2)*pi/180) + point.getZ();

            rCoord[in + 9] = r*cos(sumAngle1*pi/180)*cos((sumAngle2+dAngle2)*pi/180) + point.getX();
            rCoord[in + 10] = r*sin(sumAngle1*pi/180) + point.getY();
            rCoord[in + 11] = r*cos(sumAngle1*pi/180)*sin((sumAngle2+dAngle2)*pi/180) + point.getZ();

            in += 12;
            sumAngle1 = sumAngle1 + dAngle1;
        }
        sumAngle2 = sumAngle2 + dAngle2;
    }
    sumAngle1 = 0;
    dAngle1 = 1/(double)numStepRad;
    dAngle2 = 1/(double)numStepVer;
    in = 0;
    for(i=0; i<numStepRad; i++)
    {
        sumAngle2 = 0;
        for(j=0; j<numStepVer; j++)
        {
            tCoord[in] = sumAngle1;
            tCoord[in + 1] = sumAngle2;

            tCoord[in + 2] = sumAngle1;
            tCoord[in + 3] = sumAngle2 + dAngle2;

            tCoord[in + 4] = sumAngle1 + dAngle1;
            tCoord[in + 5] = sumAngle2 + dAngle2;

            tCoord[in + 6] = sumAngle1 + dAngle1;
            tCoord[in + 7] = sumAngle2;
            in += 8;
            sumAngle2 = sumAngle2 + dAngle2;
        }
        sumAngle1 = sumAngle1 + dAngle1;
        if(sumAngle1 > 1) sumAngle1 = 1;
    }
    for(i=0; i<nRect*12; i+=3)
    {
        nCoord[i] = (rCoord[i] - point.getX())/r;
        nCoord[i + 1] = (rCoord[i + 1] - point.getY())/r;
        nCoord[i + 2] = (rCoord[i + 2] - point.getZ())/r;
    }
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, rCoord);
    glTexCoordPointer(2, GL_DOUBLE, 0, tCoord);
    glNormalPointer(GL_DOUBLE, 0, nCoord);
    glDrawArrays(GL_QUADS, 0, nRect*4);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    delete[] rCoord;
    delete[] tCoord;
    delete[] nCoord;
}

void GetConusVertex(GLdouble r, MathPoint3D begin, MathPoint3D end, std::vector<double> &points, std::vector<double> &normals, unsigned int numSurfaces)
{
    uint32_t i;
    MathVector3D vector(begin, end), vector1, vector2, normal;
    MathPoint3D point1, point2;

    for(i=0; i<numSurfaces; i++)
    {
        points.push_back(end.getX());
        points.push_back(end.getY());
        points.push_back(end.getZ());
        point1 = GetPointFromCircle(r, vector, i*360.0/numSurfaces);
        points.push_back(point1.getX());
        points.push_back(point1.getY());
        points.push_back(point1.getZ());
        point2 = GetPointFromCircle(r, vector, (i + 1)*360.0/numSurfaces);
        points.push_back(point2.getX());
        points.push_back(point2.getY());
        points.push_back(point2.getZ());
        vector1 = MathVector3D(end, point1);
        vector2 = MathVector3D(end, point2);
        normal = (vector1^vector2).unitVector();
        normals.push_back(normal.getX());
        normals.push_back(normal.getY());
        normals.push_back(normal.getZ());
        normals.push_back(normal.getX());
        normals.push_back(normal.getY());
        normals.push_back(normal.getZ());
        normals.push_back(normal.getX());
        normals.push_back(normal.getY());
        normals.push_back(normal.getZ());
    }
}


//EOF
