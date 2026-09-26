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
#include <QFutureWatcher>
#include <atomic>
#include <vector>

/* Result of the background mesh(+solve) pipeline handed back to the GUI thread. */
struct PipelineResult
{
    bool                ok = false;
    QString             err;
    size_t              nodes = 0;
    size_t              tets = 0;
    bool                solved = false;   /* true if the eigen-solve ran */
    bool                cancelled = false;/* aborted by the user            */
    std::vector<double> k2;               /* eigenvalues, when solved */
};

/* Parametric resonator GUI: choose cavity (cylinder/box) and core (none/cylinder/
 * box/spiral with elliptical cross-section), mesh with Gmsh/OpenCASCADE, solve the
 * edge-element eigenproblem, and list the resonant frequencies. A VTK view shows
 * three toggleable layers: solids, fields, mesh. Meshing and solving run on a
 * background thread so the GUI stays responsive. */
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
    QCheckBox   *solidsCheck, *fieldsCheck, *meshCheck, *mesh3dCheck;
    QPushButton *buildMeshButton, *computeButton, *abortButton, *saveImageButton;
    QListWidget *resultsList;
    QTextEdit   *console;

    ResonatorSpec readSpec();
    void          log(const QString& s);
    void          showModeField(int internalIndex);
    void          setBusy(bool busy);
    void          finishAutoRun();   /* screenshot + print + quit, after async compute */

    /* Background pipeline (runs off the GUI thread). */
    QFutureWatcher<PipelineResult> *worker;
    bool                            solveAfterMesh;   /* current job: mesh-only vs mesh+solve */
    bool                            autoRunPending;   /* finish the headless self-test on completion */
    std::atomic<bool>               abortReq;         /* raised by Abort, polled by the worker/engine */

private slots:
    void buildMeshSlot();
    void computeSlot();
    void abortSlot();        /* request cancellation of the running job */
    void workerFinished();   /* GUI-thread handler for pipeline results */
    void cavityChanged();
    void coreChanged();
    void layerToggled();
    void modeSelected(int row);
    void saveImageSlot();
    void autoRun();   /* headless self-test: compute, screenshot, print, quit */
};

#endif // MAINDIALOG_H
