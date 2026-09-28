#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>
#include "MainDialog.h"
#include "MathEighValues.h"

int main(int argc, char** argv)
{
    EigenSolverRuntime solverRuntime(&argc, &argv);   /* MPI + hypre, for the AME solver */
    /* Must be set before the QApplication for QVTKOpenGLNativeWidget. */
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    MainDialog dialog;
    dialog.show();
    return app.exec();
}
