#ifndef MAINDIALOG_H
#define MAINDIALOG_H

#include "MicroEngine.h"
#include "MathGrapher.h"
#include "CsgGmshMesher.h"
#include <QDialog>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QListWidget>
#include <QTextEdit>
#include <QGroupBox>

/* Parametric resonator GUI: choose cavity (cylinder/box) and core (none/cylinder/
 * box/spiral with elliptical cross-section), mesh with Gmsh/OpenCASCADE, solve the
 * edge-element eigenproblem, and list the resonant frequencies. The mesh is shown
 * in the OpenGL view. */
class MainDialog : public QDialog
{
    Q_OBJECT
public:
    MainDialog();

private:
    MicroEngine *engine;
    MathGrapher *grapher;

    QComboBox   *cavityCombo, *coreCombo;
    /* cavity */
    QLineEdit   *cavRadiusEdit, *cavHeightEdit;      /* cylinder */
    QLineEdit   *cavAEdit, *cavBEdit, *cavDEdit;     /* box      */
    /* core */
    QLineEdit   *coreRadiusEdit, *coreHeightEdit;    /* cylinder */
    QLineEdit   *coreAEdit, *coreBEdit, *coreDEdit;  /* box      */
    QLineEdit   *helixREdit, *pitchEdit, *turnsEdit, *ellAEdit, *ellBEdit; /* spiral */
    QGroupBox   *cavCylBox, *cavBoxBox, *coreCylBox, *coreBoxBox, *coreSpiralBox;

    QLineEdit   *meshSizeEdit, *numModesEdit, *penaltyEdit;
    QPushButton *computeButton, *saveImageButton;
    QListWidget *resultsList;
    QTextEdit   *console;

    ResonatorSpec readSpec();
    void          log(const QString& s);

private slots:
    void computeSlot();
    void cavityChanged();
    void coreChanged();
    void saveImageSlot();
    void autoRun();   /* headless self-test: compute, screenshot, print, quit */
};

#endif // MAINDIALOG_H
