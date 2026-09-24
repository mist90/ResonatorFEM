#ifndef GLFUNCTIONS_H
#define GLFUNCTIONS_H
#include <GL/gl.h>
#include <vector>
#include "MathPoint3D.h"
#include "MathGeomObjects.h"


void GetSphereVertex(GLdouble r, MathPoint3D point, std::vector<double>& points, std::vector<double> &normals, unsigned int numStepVer=15, unsigned int numStepRad=15);
void DrawSphereTexture(GLdouble r, MathPoint3D point, GLuint numTextures, GLdouble xozAngle=0, unsigned int numStepVer=30,unsigned int numStepRad=30);

void GetConusVertex(GLdouble r, MathPoint3D begin, MathPoint3D end, std::vector<double>& points, std::vector<double> &normals, unsigned int numSurfaces=15);

#endif // GLFUNCTIONS_H
