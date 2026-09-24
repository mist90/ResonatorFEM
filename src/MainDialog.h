#ifndef MAINDIALOG_H
#define MAINDIALOG_H
#include "MicroEngine.h"
#include "MathGrapher.h"
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QDialog>
#include <QDoubleValidator>
#include <QLabel>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextEdit>

#define ACTION_START_GENERATE           (0)
#define ACTION_START_OBJECT_GENERATE    (1)
#define ACTION_START_PORT_GENERATE      (2)
#define ACTION_START_DRAW_OBJECTS       (3)
#define ACTION_END_GENERATE             (4)
#define ACTION_START_CALCULATE          (5)
#define ACTION_START_CREATE_MATRIX      (6)
#define ACTION_START_SOLVE_MATRIX       (7)
#define ACTION_END_CALCULATE            (8)
#define ACTION_END_CALCULATE_ERROR      (9)
#define ACTION_END_DRAW_VECTORS         (10)

class MainDialog: public QDialog
{
    Q_OBJECT
public:
    MainDialog();
private:
    MicroEngine *engine;
    MathGrapher *grapher;
    QPushButton *genButton, *calcButton, *drawVectorButton;
    QPushButton *setModelButton, *setVectorButton, *setAllButton, *switchTetraedrButton;
    QPushButton *saveImageButton;
    QLineEdit   *numEighValueEdit;
    QLineEdit   *phaseEdit;
    QLineEdit   *stepDrawEdit;
    QLineEdit   *lenWaveEdit, *freqEdit;
    QLineEdit   *lenStructEdit, *radCilindrEdit, *radSpiralEdit, *radDepthSpiralEdit, *numPeriodEdit;
    QTextEdit   *consolOut;
    bool        visibleTetraedr;
    void addTextConsol(uint32_t numAction);
    void panelEnable(bool enable);
private slots:
    void generateSlot();
    void generateObjectsSlot();
    void generatePortsSlot();
    void drawObjects();
    void postGenerateSlot();
    void calculateSlot();
    void createMatrixSlot();
    void solveMatrixSlot();
    void postCalculateSlot(bool ok);
    void drawVectorSlot();
    void setModelSlot();
    void setVectorSlot();
    void setAllSlot();
    void switchTetraedrSlot();
    void saveImageSlot();
};

#endif // MAINWIDGET_H
