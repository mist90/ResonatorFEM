#include <QApplication>
#include <QTextCodec>
#include <QtGlobal>
#include "MainDialog.h"


int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    MainDialog *dialog;

#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
#else
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
    QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8"));
#endif
    QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8"));
    dialog = new MainDialog();
    dialog->show();
    return app.exec();
}



//EOF
