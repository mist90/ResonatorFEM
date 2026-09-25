#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>
#include "MainDialog.h"

int main(int argc, char** argv)
{
    /* Must be set before the QApplication for QVTKOpenGLNativeWidget. */
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    MainDialog dialog;
    dialog.show();
    return app.exec();
}
