#ifndef MAINDIALOG_H
#define MAINDIALOG_H

#include "MicroEngine.h"
#include "ResonatorView.h"
#include "CsgGmshMesher.h"
#include <QDialog>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QListWidget>
#include <QTextEdit>
#include <QGroupBox>

/* Parametric resonator GUI: choose cavity (cylinder/box) and core (none/cylinder/
 * box/spiral with elliptical cross-section), mesh with Gmsh/OpenCASCADE, solve the
 * edge-element eigenproblem, and list the resonant frequencies. A VTK view shows
 * three toggleable layers: solids, fields, mesh. */
class MainDialog : public QDialog
{
    Q_OBJECT
public:
    MainDialog();

private:
    MicroEngine   *engine;
    ResonatorView *view;
    FemMesh        mesh;

    QComboBox   *cavityCombo, *coreCombo;
    QLineEdit   *cavRadiusEdit, *cavHeightEdit;
    QLineEdit   *cavAEdit, *cavBEdit, *cavDEdit;
    QLineEdit   *coreRadiusEdit, *coreHeightEdit;
    QLineEdit   *coreAEdit, *coreBEdit, *coreDEdit;
    QLineEdit   *helixREdit, *pitchEdit, *turnsEdit, *ellAEdit, *ellBEdit;
    QGroupBox   *cavCylBox, *cavBoxBox, *coreCylBox, *coreBoxBox, *coreSpiralBox;

    QLineEdit   *meshSizeEdit, *numModesEdit, *penaltyEdit;
    QCheckBox   *solidsCheck, *fieldsCheck, *meshCheck;
    QPushButton *computeButton, *saveImageButton;
    QListWidget *resultsList;
    QTextEdit   *console;

    ResonatorSpec readSpec();
    void          log(const QString& s);
    void          showModeField(int internalIndex);

private slots:
    void computeSlot();
    void cavityChanged();
    void coreChanged();
    void layerToggled();
    void modeSelected(int row);
    void saveImageSlot();
    void autoRun();   /* headless self-test: compute, screenshot, print, quit */
};

#endif // MAINDIALOG_H
