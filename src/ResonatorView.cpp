#include "ResonatorView.h"

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkTetra.h>
#include <vtkCellType.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkLookupTable.h>
#include <cmath>
#include <vector>
#include <algorithm>

ResonatorView::ResonatorView(QWidget* parent)
    : QVTKOpenGLNativeWidget(parent), modelDiagonal(1.0)
{
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renWin =
        vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    setRenderWindow(renWin);

    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(1.0, 1.0, 1.0);
    renderer->SetBackground2(0.88, 0.90, 0.96);
    renderer->GradientBackgroundOn();
    renWin->AddRenderer(renderer);

    scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetTitle("|E| (a.u.)");
    scalarBar->SetNumberOfLabels(4);
    scalarBar->SetVisibility(false);
    renderer->AddActor2D(scalarBar);
}

void ResonatorView::setMesh(const FemMesh& mesh)
{
    if (solidsActor) renderer->RemoveActor(solidsActor);
    if (meshActor)   renderer->RemoveActor(meshActor);
    if (mesh3dActor) renderer->RemoveActor(mesh3dActor);
    clearField();

    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(mesh.nodes.size());
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
    for (vtkIdType i = 0; i < (vtkIdType)mesh.nodes.size(); ++i) {
        pts->SetPoint(i, mesh.nodes[i][0], mesh.nodes[i][1], mesh.nodes[i][2]);
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], mesh.nodes[i][k]);
            hi[k] = std::max(hi[k], mesh.nodes[i][k]);
        }
    }
    modelDiagonal = std::sqrt((hi[0]-lo[0])*(hi[0]-lo[0]) + (hi[1]-lo[1])*(hi[1]-lo[1]) +
                              (hi[2]-lo[2])*(hi[2]-lo[2]));
    if (modelDiagonal <= 0.0) modelDiagonal = 1.0;

    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->SetPoints(pts);
    ug->Allocate(mesh.tets.size());
    for (const auto& t : mesh.tets) {
        vtkIdType ids[4] = { (vtkIdType)t[0], (vtkIdType)t[1], (vtkIdType)t[2], (vtkIdType)t[3] };
        ug->InsertNextCell(VTK_TETRA, 4, ids);
    }

    /* Solids: exterior surface, semi-transparent. */
    vtkSmartPointer<vtkGeometryFilter> geo = vtkSmartPointer<vtkGeometryFilter>::New();
    geo->SetInputData(ug);
    vtkSmartPointer<vtkPolyDataMapper> smap = vtkSmartPointer<vtkPolyDataMapper>::New();
    smap->SetInputConnection(geo->GetOutputPort());
    smap->ScalarVisibilityOff();
    solidsActor = vtkSmartPointer<vtkActor>::New();
    solidsActor->SetMapper(smap);
    solidsActor->GetProperty()->SetColor(0.78, 0.80, 0.86);
    solidsActor->GetProperty()->SetOpacity(0.30);
    renderer->AddActor(solidsActor);

    /* Mesh: tetrahedral wireframe. */
    vtkSmartPointer<vtkDataSetMapper> mmap = vtkSmartPointer<vtkDataSetMapper>::New();
    mmap->SetInputData(ug);
    mmap->ScalarVisibilityOff();
    meshActor = vtkSmartPointer<vtkActor>::New();
    meshActor->SetMapper(mmap);
    meshActor->GetProperty()->SetRepresentationToWireframe();
    meshActor->GetProperty()->SetColor(0.15, 0.35, 0.15);
    meshActor->GetProperty()->SetLineWidth(0.5);
    meshActor->SetVisibility(false);   /* mesh layer off by default */
    renderer->AddActor(meshActor);

    /* 3D mesh: shaded exterior surface with visible cell edges (solid look). */
    vtkSmartPointer<vtkDataSetMapper> m3 = vtkSmartPointer<vtkDataSetMapper>::New();
    m3->SetInputData(ug);
    m3->ScalarVisibilityOff();
    mesh3dActor = vtkSmartPointer<vtkActor>::New();
    mesh3dActor->SetMapper(m3);
    mesh3dActor->GetProperty()->SetRepresentationToSurface();
    mesh3dActor->GetProperty()->EdgeVisibilityOn();
    mesh3dActor->GetProperty()->SetColor(0.80, 0.82, 0.86);
    mesh3dActor->GetProperty()->SetEdgeColor(0.20, 0.30, 0.20);
    mesh3dActor->SetVisibility(false);   /* off by default */
    renderer->AddActor(mesh3dActor);

    resetView();
}

void ResonatorView::setField(const std::vector<std::array<double, 3> >& points,
                             const std::vector<std::array<double, 3> >& vectors)
{
    clearField();
    if (points.empty() || points.size() != vectors.size()) return;

    /* Magnitudes. A few tetrahedra near curved walls can be slivers, where the
     * Whitney-basis field blows up (~1/volume). Use a robust reference scale
     * (95th percentile) for colour and arrow size, and drop the wild outliers,
     * so those bad elements don't dwarf the physical field. */
    std::vector<double> mags(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        mags[i] = std::sqrt(vectors[i][0]*vectors[i][0] + vectors[i][1]*vectors[i][1] +
                            vectors[i][2]*vectors[i][2]);
    std::vector<double> sorted(mags);
    std::sort(sorted.begin(), sorted.end());
    double robustMax = sorted[(size_t)(0.95 * (sorted.size() - 1))];
    if (robustMax <= 0.0) robustMax = sorted.back();
    if (robustMax <= 0.0) return;
    const double outlierCut = 6.0 * robustMax;   /* drop clearly non-physical samples */

    vtkSmartPointer<vtkPoints>      p   = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> vec = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> mag = vtkSmartPointer<vtkDoubleArray>::New();
    vec->SetNumberOfComponents(3); vec->SetName("field");
    mag->SetNumberOfComponents(1); mag->SetName("mag");
    for (size_t i = 0; i < points.size(); ++i) {
        if (mags[i] > outlierCut) continue;
        p->InsertNextPoint(points[i][0], points[i][1], points[i][2]);
        vec->InsertNextTuple3(vectors[i][0], vectors[i][1], vectors[i][2]);
        mag->InsertNextValue(mags[i]);
    }

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(p);
    pd->GetPointData()->SetVectors(vec);
    pd->GetPointData()->SetScalars(mag);

    vtkSmartPointer<vtkArrowSource> arrow = vtkSmartPointer<vtkArrowSource>::New();
    vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
    glyph->SetInputData(pd);
    glyph->SetSourceConnection(arrow->GetOutputPort());
    glyph->SetVectorModeToUseVector();
    glyph->SetScaleModeToScaleByVector();
    glyph->SetColorModeToColorByScalar();
    /* Clamp the scaling to [0, robustMax] so arrows above the 95th percentile
     * saturate at one size instead of dominating. */
    glyph->ClampingOn();
    glyph->SetRange(0.0, robustMax);
    glyph->SetScaleFactor(0.06 * modelDiagonal);
    glyph->OrientOn();

    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetHueRange(0.667, 0.0);   /* blue (low) -> red (high) */
    lut->SetTableRange(0.0, robustMax);
    lut->Build();

    vtkSmartPointer<vtkPolyDataMapper> map = vtkSmartPointer<vtkPolyDataMapper>::New();
    map->SetInputConnection(glyph->GetOutputPort());
    map->SetLookupTable(lut);
    map->SetScalarRange(0.0, robustMax);

    fieldsActor = vtkSmartPointer<vtkActor>::New();
    fieldsActor->SetMapper(map);
    renderer->AddActor(fieldsActor);

    scalarBar->SetLookupTable(lut);
    scalarBar->SetVisibility(true);
    renderWindow()->Render();
}

void ResonatorView::clearField()
{
    if (fieldsActor) { renderer->RemoveActor(fieldsActor); fieldsActor = nullptr; }
    if (scalarBar)   scalarBar->SetVisibility(false);
}

void ResonatorView::setLayerVisible(Layer layer, bool visible)
{
    switch (layer) {
    case SOLIDS: if (solidsActor) solidsActor->SetVisibility(visible); break;
    case MESH:   if (meshActor)   meshActor->SetVisibility(visible);   break;
    case MESH3D: if (mesh3dActor) mesh3dActor->SetVisibility(visible); break;
    case FIELDS:
        if (fieldsActor) fieldsActor->SetVisibility(visible);
        scalarBar->SetVisibility(visible && fieldsActor != nullptr);
        break;
    }
    if (renderWindow()) renderWindow()->Render();
}

void ResonatorView::resetView()
{
    renderer->ResetCamera();
    if (renderWindow()) renderWindow()->Render();
}
