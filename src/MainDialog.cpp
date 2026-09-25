#include "MainDialog.h"
#include <QApplication>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QLabel>
#include <QDoubleValidator>
#include <QIntValidator>
#include <QLocale>
#include <QFileDialog>
#include <QListWidgetItem>
#include <QTimer>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

static const double C_LIGHT = 299792458.0;   /* m/s */

static QLineEdit* field(const QString& val)
{
    QLineEdit* e = new QLineEdit(val);
    QDoubleValidator* v = new QDoubleValidator(0.0, 1e9, 6, e);
    /* Force '.' as the decimal separator regardless of system locale, and plain
     * (non-scientific) notation. Otherwise on locales that use ',' the validator
     * rejects '.' and the fields feel un-editable. */
    v->setNotation(QDoubleValidator::StandardNotation);
    v->setLocale(QLocale::c());
    e->setValidator(v);
    return e;
}

MainDialog::MainDialog()
{
    setWindowTitle("ResonatorFEM");
    view   = new ResonatorView();
    engine = new MicroEngine();

    cavityCombo = new QComboBox();
    cavityCombo->addItems({ "Cylinder", "Rectangular (box)" });
    coreCombo = new QComboBox();
    coreCombo->addItems({ "None", "Cylinder", "Rectangular (box)", "Spiral (elliptical)" });

    cavRadiusEdit = field("1.0"); cavHeightEdit = field("2.0");
    cavCylBox = new QGroupBox("Cylinder cavity");
    { QFormLayout* f = new QFormLayout(cavCylBox);
      f->addRow("radius (m)", cavRadiusEdit); f->addRow("height (m)", cavHeightEdit); }

    cavAEdit = field("1.0"); cavBEdit = field("0.8"); cavDEdit = field("1.5");
    cavBoxBox = new QGroupBox("Box cavity");
    { QFormLayout* f = new QFormLayout(cavBoxBox);
      f->addRow("a — x (m)", cavAEdit); f->addRow("b — y (m)", cavBEdit); f->addRow("d — z (m)", cavDEdit); }

    coreRadiusEdit = field("0.3"); coreHeightEdit = field("2.0");
    coreCylBox = new QGroupBox("Cylinder core");
    { QFormLayout* f = new QFormLayout(coreCylBox);
      f->addRow("radius (m)", coreRadiusEdit); f->addRow("height (m)", coreHeightEdit); }

    coreAEdit = field("0.3"); coreBEdit = field("0.3"); coreDEdit = field("1.0");
    coreBoxBox = new QGroupBox("Box core");
    { QFormLayout* f = new QFormLayout(coreBoxBox);
      f->addRow("a — x (m)", coreAEdit); f->addRow("b — y (m)", coreBEdit); f->addRow("d — z (m)", coreDEdit); }

    helixREdit = field("0.5"); pitchEdit = field("0.5"); turnsEdit = field("2.0");
    ellAEdit = field("0.12");  ellBEdit = field("0.12");
    coreSpiralBox = new QGroupBox("Spiral core (elliptical cross-section)");
    { QFormLayout* f = new QFormLayout(coreSpiralBox);
      f->addRow("helix radius (m)", helixREdit); f->addRow("pitch/turn (m)", pitchEdit);
      f->addRow("turns (–)", turnsEdit); f->addRow("semi-axis radial (m)", ellAEdit);
      f->addRow("semi-axis axial (m)", ellBEdit); }

    meshSizeEdit = field("0.2");
    numModesEdit = new QLineEdit("8"); numModesEdit->setValidator(new QIntValidator(1, 100));
    penaltyEdit  = field("50");

    solidsCheck = new QCheckBox("Solids"); solidsCheck->setChecked(true);
    fieldsCheck = new QCheckBox("Fields"); fieldsCheck->setChecked(true);
    meshCheck   = new QCheckBox("Mesh");   meshCheck->setChecked(false);

    computeButton   = new QPushButton("Compute modes");
    saveImageButton = new QPushButton("Save image");
    resultsList = new QListWidget();
    console = new QTextEdit(); console->setReadOnly(true); console->setMaximumHeight(110);

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
    solveForm->addRow("mesh size (m)", meshSizeEdit);
    solveForm->addRow("# modes", numModesEdit);
    solveForm->addRow("penalty factor", penaltyEdit);
    controls->addLayout(solveForm);
    controls->addWidget(computeButton);
    QLabel* unitsNote = new QLabel(
        "Lengths in metres · frequencies in MHz · field |E| in\n"
        "arbitrary units (an eigenmode is defined up to a scale).");
    unitsNote->setStyleSheet("color: #555; font-size: 10px;");
    controls->addWidget(unitsNote);
    controls->addWidget(new QLabel("Resonant frequencies (MHz):"));
    controls->addWidget(resultsList);
    controls->addWidget(console);

    QHBoxLayout* layerRow = new QHBoxLayout();
    layerRow->addWidget(new QLabel("Layers:"));
    layerRow->addWidget(solidsCheck);
    layerRow->addWidget(fieldsCheck);
    layerRow->addWidget(meshCheck);
    layerRow->addStretch(1);
    layerRow->addWidget(saveImageButton);

    QVBoxLayout* viewCol = new QVBoxLayout();
    viewCol->addWidget(view, 1);
    viewCol->addLayout(layerRow);

    QHBoxLayout* root = new QHBoxLayout(this);
    QWidget* controlsW = new QWidget(); controlsW->setLayout(controls);
    controlsW->setMaximumWidth(340);
    root->addWidget(controlsW);
    root->addLayout(viewCol, 1);

    connect(cavityCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(cavityChanged()));
    connect(coreCombo,   SIGNAL(currentIndexChanged(int)), this, SLOT(coreChanged()));
    connect(computeButton,   SIGNAL(clicked()), this, SLOT(computeSlot()));
    connect(saveImageButton, SIGNAL(clicked()), this, SLOT(saveImageSlot()));
    connect(solidsCheck, SIGNAL(toggled(bool)), this, SLOT(layerToggled()));
    connect(fieldsCheck, SIGNAL(toggled(bool)), this, SLOT(layerToggled()));
    connect(meshCheck,   SIGNAL(toggled(bool)), this, SLOT(layerToggled()));
    connect(resultsList, SIGNAL(currentRowChanged(int)), this, SLOT(modeSelected(int)));

    cavityChanged();
    coreChanged();
    resize(1060, 660);

    if (std::getenv("RESONATOR_AUTORUN")) QTimer::singleShot(400, this, SLOT(autoRun()));
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

void MainDialog::layerToggled()
{
    view->setLayerVisible(ResonatorView::SOLIDS, solidsCheck->isChecked());
    view->setLayerVisible(ResonatorView::FIELDS, fieldsCheck->isChecked());
    view->setLayerVisible(ResonatorView::MESH,   meshCheck->isChecked());
}

void MainDialog::showModeField(int internalIndex)
{
    if (internalIndex < 0) return;
    engine->calculateBasisKoef((uint32_t)internalIndex);
    std::vector<std::array<double, 3> > pts, vecs;
    engine->sampleFieldAtCentroids(pts, vecs);
    view->setField(pts, vecs);
    view->setLayerVisible(ResonatorView::FIELDS, fieldsCheck->isChecked());
}

void MainDialog::modeSelected(int row)
{
    if (row < 0 || row >= resultsList->count()) return;
    showModeField(resultsList->item(row)->data(Qt::UserRole).toInt());
}

void MainDialog::computeSlot()
{
    computeButton->setEnabled(false);
    resultsList->clear();
    ResonatorSpec spec = readSpec();

    log("Meshing (Gmsh/OpenCASCADE)...");
    std::string err;
    if (!CsgGmshMesher::buildResonator(spec, mesh, &err)) {
        log("Mesh FAILED: " + QString::fromStdString(err));
        computeButton->setEnabled(true);
        return;
    }
    log(QString("Mesh: %1 nodes, %2 tets").arg(mesh.nodes.size()).arg(mesh.tets.size()));
    view->setMesh(mesh);
    layerToggled();

    engine->setResonatorMode(true);
    engine->setGradDivPenalty(penaltyEdit->text().toDouble());
    int nModes = numModesEdit->text().toInt();
    engine->setResonatorSolveTarget(0.001, nModes > 0 ? nModes : 8);
    if (!engine->generateFromMesh(mesh)) {
        log("generateFromMesh failed");
        computeButton->setEnabled(true);
        return;
    }
    log("Solving eigenproblem...");
    if (!engine->calculateSync()) {
        log("Solve FAILED");
        computeButton->setEnabled(true);
        return;
    }

    std::vector<double> k2;
    engine->getEighValues(k2);
    std::vector<std::pair<double, int> > modes;
    for (int i = 0; i < (int)k2.size(); ++i) modes.push_back(std::make_pair(k2[i], i));
    std::sort(modes.begin(), modes.end());
    for (const auto& m : modes) {
        if (m.first <= 0.0) continue;
        double f = C_LIGHT * std::sqrt(m.first) / (2.0 * M_PI) / 1.0e6;
        QListWidgetItem* it = new QListWidgetItem(
            QString("%1 MHz   (k^2 = %2)").arg(f, 0, 'f', 3).arg(m.first, 0, 'f', 5));
        it->setData(Qt::UserRole, m.second);   /* internal mode index for basis reconstruction */
        resultsList->addItem(it);
    }
    if (resultsList->count() > 0) resultsList->setCurrentRow(0);   /* show fundamental field */
    log("Done.");
    computeButton->setEnabled(true);
}

void MainDialog::saveImageSlot()
{
    QString fn = QFileDialog::getSaveFileName(this, "Save image", "resonator.png", "PNG (*.png)");
    if (!fn.isEmpty()) view->grabFramebuffer().save(fn);
}

void MainDialog::autoRun()
{
    if (const char* c = std::getenv("RESONATOR_CAV"))  cavityCombo->setCurrentIndex(std::atoi(c));
    if (const char* c = std::getenv("RESONATOR_CORE")) coreCombo->setCurrentIndex(std::atoi(c));
    computeSlot();
    if (const char* shot = std::getenv("RESONATOR_SHOT")) {
        view->grabFramebuffer().save(QString(shot));
        this->grab().save(QString(shot) + ".ui.png");   /* full window (controls) */
    }
    std::printf("=== AUTORUN results (%d modes) ===\n", resultsList->count());
    for (int i = 0; i < resultsList->count(); ++i)
        std::printf("  %s\n", resultsList->item(i)->text().toStdString().c_str());
    std::fflush(stdout);
    qApp->quit();
}
