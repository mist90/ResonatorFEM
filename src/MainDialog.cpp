#include "MainDialog.h"
#include <QApplication>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QLabel>
#include <QDoubleValidator>
#include <QIntValidator>
#include <QFileDialog>
#include <QTimer>
#include <QPixmap>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <vector>

static const double C_LIGHT = 299792458.0;   /* m/s */

static QLineEdit* field(const QString& val)
{
    QLineEdit* e = new QLineEdit(val);
    e->setValidator(new QDoubleValidator(0.0, 1e9, 6));
    return e;
}

MainDialog::MainDialog()
{
    setWindowTitle("ResonatorFEM");
    grapher = new MathGrapher();
    engine  = new MicroEngine(*grapher);

    /* ---- shape selectors ---- */
    cavityCombo = new QComboBox();
    cavityCombo->addItems({ "Cylinder", "Rectangular (box)" });
    coreCombo = new QComboBox();
    coreCombo->addItems({ "None", "Cylinder", "Rectangular (box)", "Spiral (elliptical)" });

    /* ---- cavity parameter groups ---- */
    cavRadiusEdit = field("1.0"); cavHeightEdit = field("2.0");
    cavCylBox = new QGroupBox("Cylinder cavity");
    { QFormLayout* f = new QFormLayout(cavCylBox);
      f->addRow("radius", cavRadiusEdit); f->addRow("height", cavHeightEdit); }

    cavAEdit = field("1.0"); cavBEdit = field("0.8"); cavDEdit = field("1.5");
    cavBoxBox = new QGroupBox("Box cavity");
    { QFormLayout* f = new QFormLayout(cavBoxBox);
      f->addRow("a (x)", cavAEdit); f->addRow("b (y)", cavBEdit); f->addRow("d (z)", cavDEdit); }

    /* ---- core parameter groups ---- */
    coreRadiusEdit = field("0.3"); coreHeightEdit = field("2.0");
    coreCylBox = new QGroupBox("Cylinder core");
    { QFormLayout* f = new QFormLayout(coreCylBox);
      f->addRow("radius", coreRadiusEdit); f->addRow("height", coreHeightEdit); }

    coreAEdit = field("0.3"); coreBEdit = field("0.3"); coreDEdit = field("1.0");
    coreBoxBox = new QGroupBox("Box core");
    { QFormLayout* f = new QFormLayout(coreBoxBox);
      f->addRow("a (x)", coreAEdit); f->addRow("b (y)", coreBEdit); f->addRow("d (z)", coreDEdit); }

    helixREdit = field("0.5"); pitchEdit = field("0.5"); turnsEdit = field("2.0");
    ellAEdit = field("0.12");  ellBEdit = field("0.12");
    coreSpiralBox = new QGroupBox("Spiral core (elliptical cross-section)");
    { QFormLayout* f = new QFormLayout(coreSpiralBox);
      f->addRow("helix radius", helixREdit); f->addRow("pitch/turn", pitchEdit);
      f->addRow("turns", turnsEdit); f->addRow("ellipse semi-axis (radial)", ellAEdit);
      f->addRow("ellipse semi-axis (axial)", ellBEdit); }

    /* ---- solve parameters ---- */
    meshSizeEdit = field("0.2");
    numModesEdit = new QLineEdit("8"); numModesEdit->setValidator(new QIntValidator(1, 100));
    penaltyEdit  = field("50");

    computeButton   = new QPushButton("Compute modes");
    saveImageButton = new QPushButton("Save image");
    resultsList = new QListWidget();
    console = new QTextEdit(); console->setReadOnly(true);
    console->setMaximumHeight(120);

    /* ---- layout: left controls, right GL view ---- */
    QVBoxLayout* controls = new QVBoxLayout();
    QFormLayout* shapeForm = new QFormLayout();
    shapeForm->addRow("Cavity", cavityCombo);
    shapeForm->addRow("Core",   coreCombo);
    controls->addLayout(shapeForm);
    controls->addWidget(cavCylBox);
    controls->addWidget(cavBoxBox);
    controls->addWidget(coreCylBox);
    controls->addWidget(coreBoxBox);
    controls->addWidget(coreSpiralBox);
    QFormLayout* solveForm = new QFormLayout();
    solveForm->addRow("mesh size", meshSizeEdit);
    solveForm->addRow("# modes", numModesEdit);
    solveForm->addRow("penalty factor", penaltyEdit);
    controls->addLayout(solveForm);
    controls->addWidget(computeButton);
    controls->addWidget(new QLabel("Resonant frequencies (MHz):"));
    controls->addWidget(resultsList);
    controls->addWidget(console);

    QVBoxLayout* viewCol = new QVBoxLayout();
    viewCol->addWidget(grapher, 1);
    viewCol->addWidget(saveImageButton);

    QHBoxLayout* root = new QHBoxLayout(this);
    QWidget* controlsW = new QWidget(); controlsW->setLayout(controls);
    controlsW->setMaximumWidth(340);
    root->addWidget(controlsW);
    root->addLayout(viewCol, 1);

    connect(cavityCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(cavityChanged()));
    connect(coreCombo,   SIGNAL(currentIndexChanged(int)), this, SLOT(coreChanged()));
    connect(computeButton,   SIGNAL(clicked()), this, SLOT(computeSlot()));
    connect(saveImageButton, SIGNAL(clicked()), this, SLOT(saveImageSlot()));

    cavityChanged();
    coreChanged();
    resize(1000, 640);

    if (std::getenv("RESONATOR_AUTORUN")) QTimer::singleShot(400, this, SLOT(autoRun()));
}

void MainDialog::autoRun()
{
    computeSlot();
    if (const char* shot = std::getenv("RESONATOR_SHOT")) {
        this->grab().save(QString(shot));
        QImage gl = grapher->grabFrameBuffer();
        gl.save(QString(shot) + ".glview.png");
    }
    std::printf("=== AUTORUN results (%d modes) ===\n", resultsList->count());
    for (int i = 0; i < resultsList->count(); ++i)
        std::printf("  %s\n", resultsList->item(i)->text().toStdString().c_str());
    std::fflush(stdout);
    qApp->quit();
}

void MainDialog::cavityChanged()
{
    bool cyl = (cavityCombo->currentIndex() == ResonatorSpec::CAVITY_CYLINDER);
    cavCylBox->setVisible(cyl);
    cavBoxBox->setVisible(!cyl);
}

void MainDialog::coreChanged()
{
    int c = coreCombo->currentIndex();
    coreCylBox->setVisible(c == ResonatorSpec::CORE_CYLINDER);
    coreBoxBox->setVisible(c == ResonatorSpec::CORE_BOX);
    coreSpiralBox->setVisible(c == ResonatorSpec::CORE_SPIRAL);
}

ResonatorSpec MainDialog::readSpec()
{
    ResonatorSpec s;
    s.cavity     = cavityCombo->currentIndex();
    s.cavRadius  = cavRadiusEdit->text().toDouble();
    s.cavHeight  = cavHeightEdit->text().toDouble();
    s.cavA = cavAEdit->text().toDouble();
    s.cavB = cavBEdit->text().toDouble();
    s.cavD = cavDEdit->text().toDouble();
    s.core       = coreCombo->currentIndex();
    s.coreRadius = coreRadiusEdit->text().toDouble();
    s.coreHeight = coreHeightEdit->text().toDouble();
    s.coreA = coreAEdit->text().toDouble();
    s.coreB = coreBEdit->text().toDouble();
    s.coreD = coreDEdit->text().toDouble();
    s.helixR = helixREdit->text().toDouble();
    s.pitch  = pitchEdit->text().toDouble();
    s.turns  = turnsEdit->text().toDouble();
    s.ellA   = ellAEdit->text().toDouble();
    s.ellB   = ellBEdit->text().toDouble();
    s.meshSize = meshSizeEdit->text().toDouble();
    return s;
}

void MainDialog::log(const QString& s)
{
    console->append(s);
    QApplication::processEvents();
}

void MainDialog::computeSlot()
{
    computeButton->setEnabled(false);
    resultsList->clear();
    ResonatorSpec spec = readSpec();

    log("Meshing (Gmsh/OpenCASCADE)...");
    FemMesh mesh;
    std::string err;
    if (!CsgGmshMesher::buildResonator(spec, mesh, &err)) {
        log("Mesh FAILED: " + QString::fromStdString(err));
        computeButton->setEnabled(true);
        return;
    }
    log(QString("Mesh: %1 nodes, %2 tets").arg(mesh.nodes.size()).arg(mesh.tets.size()));

    /* View bounds from the mesh. */
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
    for (const auto& p : mesh.nodes)
        for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], p[k]); hi[k] = std::max(hi[k], p[k]); }
    grapher->clearPoints(); grapher->clearLines(); grapher->clearTetraedrs();
    grapher->setBeginPoint(MathPoint3D(lo[0], lo[1], lo[2]));
    grapher->setSizeRegion(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);

    engine->setResonatorMode(true);
    engine->setGradDivPenalty(penaltyEdit->text().toDouble());
    int nModes = numModesEdit->text().toInt();
    engine->setResonatorSolveTarget(0.001, nModes > 0 ? nModes : 8);
    if (!engine->generateFromMesh(mesh)) {
        log("generateFromMesh failed");
        computeButton->setEnabled(true);
        return;
    }
    engine->drawAllObject();     /* store mesh points + tetrahedra in the view */
    grapher->setEnable(POINTS_TYPE);
    grapher->setEnable(TETRAEDRS_TYPE);
    grapher->setEnable(DRAW);    /* enables drawing + triggers a repaint */

    log("Solving eigenproblem...");
    if (!engine->calculateSync()) {
        log("Solve FAILED");
        computeButton->setEnabled(true);
        return;
    }
    std::vector<double> k2;
    engine->getEighValues(k2);
    std::sort(k2.begin(), k2.end());
    for (double v : k2) {
        if (v <= 0.0) continue;
        double f = C_LIGHT * std::sqrt(v) / (2.0 * M_PI) / 1.0e6;
        resultsList->addItem(QString("%1  MHz   (k^2 = %2)")
                             .arg(f, 0, 'f', 3).arg(v, 0, 'f', 5));
    }
    log("Done.");
    computeButton->setEnabled(true);
}

void MainDialog::saveImageSlot()
{
    QString fn = QFileDialog::getSaveFileName(this, "Save image", "resonator.png", "PNG (*.png)");
    if (!fn.isEmpty()) grapher->saveImage(fn);
}
