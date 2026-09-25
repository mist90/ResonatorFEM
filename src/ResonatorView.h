#ifndef RESONATORVIEW_H
#define RESONATORVIEW_H

#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <array>
#include <vector>
#include "FemMesh.h"

/* VTK-based 3D view with three independently toggleable layers:
 *   SOLIDS  - the geometry: exterior (PEC) surface of the meshed domain,
 *             semi-transparent so the interior is visible
 *   MESH    - the tetrahedral mesh as a wireframe
 *   FIELDS  - the computed mode field as arrow glyphs coloured by magnitude
 */
class ResonatorView : public QVTKOpenGLNativeWidget
{
    Q_OBJECT
public:
    enum Layer { SOLIDS = 0, MESH = 1, FIELDS = 2 };

    explicit ResonatorView(QWidget* parent = nullptr);

    /* Build the solids + mesh layers from a tetrahedral mesh. Clears any field. */
    void setMesh(const FemMesh& mesh);
    /* Build the fields layer from sampled (point, vector) pairs. */
    void setField(const std::vector<std::array<double, 3> >& points,
                  const std::vector<std::array<double, 3> >& vectors);
    void clearField();

    void setLayerVisible(Layer layer, bool visible);
    void resetView();

private:
    vtkSmartPointer<vtkRenderer>        renderer;
    vtkSmartPointer<vtkActor>           solidsActor;
    vtkSmartPointer<vtkActor>           meshActor;
    vtkSmartPointer<vtkActor>           fieldsActor;
    vtkSmartPointer<vtkScalarBarActor>  scalarBar;
    double modelDiagonal;   /* for glyph scaling */
};

#endif // RESONATORVIEW_H
