#include "MainDialog.h"
#include <QApplication>
#include <float.h>

#define LBV


MainDialog::MainDialog()
{
    grapher = new MathGrapher();
    engine = new MicroEngine(*grapher);
    genButton = new QPushButton(QString("Сгенерировать сетку"));

    calcButton = new QPushButton(QString("Рассчитать поле"));
    drawVectorButton = new QPushButton(QString("Отрисовать векторы"));

    numEighValueEdit = new QLineEdit();
    QIntValidator intValidator;
    intValidator.setBottom(0);
    numEighValueEdit->setValidator(&intValidator);
    numEighValueEdit->setText(QString("0"));

    phaseEdit = new QLineEdit();
    phaseEdit->setValidator(new QDoubleValidator());
    phaseEdit->setText(QString("0"));

    lenWaveEdit = new QLineEdit();
    lenWaveEdit->setReadOnly(true);
    lenWaveEdit->setText(QString("0"));

    freqEdit = new QLineEdit();
    freqEdit->setReadOnly(true);
    freqEdit->setText(QString("0"));

    stepDrawEdit = new QLineEdit();
    QDoubleValidator validator;
    validator.setBottom(0.0);
    stepDrawEdit->setValidator(&validator);
    stepDrawEdit->setText(QString("0"));

    /* QLineEdit'ы, задающие габариты модели */
    lenStructEdit = new QLineEdit();
    lenStructEdit->setValidator(&validator);
    lenStructEdit->setText(QString("%1").arg(7.2));
    radCilindrEdit = new QLineEdit();
    radCilindrEdit->setValidator(&validator);
    radCilindrEdit->setText(QString("%1").arg(2.5));
    radSpiralEdit = new QLineEdit();
    radSpiralEdit->setValidator(&validator);
    radSpiralEdit->setText(QString("%1").arg(1.1));
    radDepthSpiralEdit = new QLineEdit();
    radDepthSpiralEdit->setValidator(&validator);
    radDepthSpiralEdit->setText(QString("%1").arg(0.2));
    numPeriodEdit = new QLineEdit();
    numPeriodEdit->setValidator(&intValidator);
    numPeriodEdit->setText(QString("%1").arg(6));

    consolOut = new QTextEdit();
    consolOut->setReadOnly(true);
    consolOut->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
    setModelButton = new QPushButton(QString("Модель"));
    setVectorButton = new QPushButton(QString("Поле"));
    setAllButton = new QPushButton(QString("Поле и модель"));
    switchTetraedrButton = new QPushButton(QString("Тетраэдры видимы: вкл"));
    visibleTetraedr = true;
    saveImageButton = new QPushButton(QString("Сохранить изображение"));

    QVBoxLayout *leftLayout = new QVBoxLayout();
    leftLayout->addWidget(grapher);

    QHBoxLayout *grapherPanelLayout = new QHBoxLayout();

    QVBoxLayout *columnModelLayout = new QVBoxLayout();
    columnModelLayout->addWidget(setModelButton);
    columnModelLayout->addWidget(switchTetraedrButton);
    grapherPanelLayout->addLayout(columnModelLayout);

    QVBoxLayout *columnVectorLayout = new QVBoxLayout();
    columnVectorLayout->addWidget(setVectorButton);

    grapherPanelLayout->addLayout(columnVectorLayout);
    grapherPanelLayout->addWidget(setAllButton);

    leftLayout->addLayout(grapherPanelLayout);
    leftLayout->addWidget(consolOut);

    QWidget *rightWidget = new QWidget();
    QVBoxLayout *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(genButton);
    rightLayout->addWidget(new QLabel(QString("Длина резонатора, мм")));
    rightLayout->addWidget(lenStructEdit);
    rightLayout->addWidget(new QLabel(QString("Радиус резонатора, мм")));
    rightLayout->addWidget(radCilindrEdit);
    rightLayout->addWidget(new QLabel(QString("Радиус спирали, мм")));
    rightLayout->addWidget(radSpiralEdit);
    rightLayout->addWidget(new QLabel(QString("Радиус сечения спирали, мм")));
    rightLayout->addWidget(radDepthSpiralEdit);
    rightLayout->addWidget(new QLabel(QString("Число периодов")));
    rightLayout->addWidget(numPeriodEdit);
    rightLayout->addSpacing(20);
    rightLayout->addWidget(calcButton);
    rightLayout->addSpacing(20);
    rightLayout->addWidget(new QLabel(QString("Номер соб. значения:")));
    rightLayout->addWidget(numEighValueEdit);
    rightLayout->addWidget(new QLabel(QString("Фаза поля, рад:")));
    rightLayout->addWidget(phaseEdit);
    rightLayout->addWidget(new QLabel(QString("Шаг отрисовки, м:")));
    rightLayout->addWidget(stepDrawEdit);
    rightLayout->addWidget(drawVectorButton);
    rightLayout->addSpacing(20);
    rightLayout->addWidget(saveImageButton);
    rightLayout->addSpacing(20);
    rightLayout->addWidget(new QLabel(QString("Длина волны, м:")));
    rightLayout->addWidget(lenWaveEdit);
    rightLayout->addWidget(new QLabel(QString("Частота, ГГц:")));
    rightLayout->addWidget(freqEdit);
    rightLayout->addStretch();
    rightWidget->setMaximumWidth(200);
    rightWidget->setLayout(rightLayout);


    QHBoxLayout *mainLayout = new QHBoxLayout();
    mainLayout->addLayout(leftLayout);
    mainLayout->addWidget(rightWidget);
    setLayout(mainLayout);

    connect(genButton, SIGNAL(clicked()), this, SLOT(generateSlot()));
    connect(engine, SIGNAL(startObjectsGenerate()), this, SLOT(generateObjectsSlot()));
    connect(engine, SIGNAL(startPortsGenerate()), this, SLOT(generatePortsSlot()));
    connect(engine, SIGNAL(startDrawObjects()), this, SLOT(drawObjects()));
    connect(engine, SIGNAL(endGenerate()), this, SLOT(postGenerateSlot()));
    connect(calcButton, SIGNAL(clicked()), this, SLOT(calculateSlot()));
    connect(engine, SIGNAL(startCreateMatrix()), this, SLOT(createMatrixSlot()));
    connect(engine, SIGNAL(startSolveMatrix()), this, SLOT(solveMatrixSlot()));
    connect(engine, SIGNAL(endCalculate(bool)), this, SLOT(postCalculateSlot(bool)));
    connect(drawVectorButton, SIGNAL(clicked()), this, SLOT(drawVectorSlot()));
    connect(setModelButton, SIGNAL(clicked()), this, SLOT(setModelSlot()));
    connect(setVectorButton, SIGNAL(clicked()), this, SLOT(setVectorSlot()));
    connect(setAllButton, SIGNAL(clicked()), this, SLOT(setAllSlot()));
    connect(switchTetraedrButton, SIGNAL(clicked()), this, SLOT(switchTetraedrSlot()));
    connect(saveImageButton, SIGNAL(clicked()), this, SLOT(saveImageSlot()));
    setModelSlot();
    setWindowFlags(Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint | Qt::WindowMinimizeButtonHint);
}

void MainDialog::panelEnable(bool enable)
{
    genButton->setEnabled(enable);
    calcButton->setEnabled(enable);
    phaseEdit->setEnabled(enable);
    stepDrawEdit->setEnabled(enable);
    drawVectorButton->setEnabled(enable);
    numEighValueEdit->setEnabled(enable);
}

void MainDialog::generateSlot()
{
    MathObject object;

    const double div = 1000;
    const double lenCilindr = lenStructEdit->text().toDouble()/div, radCilindr = radCilindrEdit->text().toDouble()/div;
    const double lenSpiral = lenStructEdit->text().toDouble()/div, radSpiral = radSpiralEdit->text().toDouble()/div;
    const double radDepthSpiral = radDepthSpiralEdit->text().toDouble()/div;
    const double stepSpiral = lenStructEdit->text().toDouble()/numPeriodEdit->text().toDouble()/div;

    MathPoint3D center(5.0, 5.0, 5.0);

    grapher->clearPoints();
    grapher->clearLines();
    grapher->clearTetraedrs();
    grapher->clearVectors();
    grapher->setRadiusVolumeLines(0.02/div);
    engine->clear();
    engine->setBeginPoint(MathPoint3D(0.0, 0.0, 0.0));
    engine->setSizeGrid(10.0, 10.0, 10.0);
    engine->setResonatorMode(true);
    addTextConsol(ACTION_START_GENERATE);

    // добавление цилиндров
    object.makeCilindr(center, radCilindr, lenCilindr);
    object.setVector(MathVector3D(1.0, 0.0, 0.0));
    engine->addObject(object, true, 0.0, 1.0);

    // добавление цилиндра для коаксиальной линии
    /*object.makeCilindr(center, radCilindr/6.0, lenCilindr, 2.0);
    object.setVector(MathVector3D(1.0, 0.0, 0.0));
    engine->addObject(object, true, SIGMA_IDEAL_METALL, 1.0);*/

    // добавление спирали для ЛБВ
    object.makeSpiral(center, radDepthSpiral, radSpiral, stepSpiral, -90.0, lenSpiral);
    object.setVector(MathVector3D(1.0, 0.0, 0.0));
    engine->addObject(object, true, SIGMA_IDEAL_METALL, 1.0);

    setEnabled(false);
    engine->genPoints();
}

void MainDialog::generateObjectsSlot()
{
    addTextConsol(ACTION_START_OBJECT_GENERATE);
}

void MainDialog::generatePortsSlot()
{
    addTextConsol(ACTION_START_PORT_GENERATE);
}

void MainDialog::drawObjects()
{
    addTextConsol(ACTION_START_DRAW_OBJECTS);
}

void MainDialog::postGenerateSlot()
{
    engine->drawAllObject();
    grapher->setEnable(DRAW);
    addTextConsol(ACTION_END_GENERATE);
    setEnabled(true);
}

void MainDialog::calculateSlot()
{
    if(!engine->isGenerate())
        QMessageBox::warning(this, QString("Ошибка"), QString("Сначала нужно сгенерировать сетку!"), QMessageBox::Ok);
    else
    {
        panelEnable(false);
        addTextConsol(ACTION_START_CALCULATE);
        engine->calculate();
    }
}

void MainDialog::createMatrixSlot()
{
    addTextConsol(ACTION_START_CREATE_MATRIX);
}

void MainDialog::solveMatrixSlot()
{
    addTextConsol(ACTION_START_SOLVE_MATRIX);
}

void MainDialog::postCalculateSlot(bool ok)
{
    if(ok) addTextConsol(ACTION_END_CALCULATE);
    else addTextConsol(ACTION_END_CALCULATE_ERROR);
    panelEnable(true);
}

void MainDialog::drawVectorSlot()
{
    uint32_t i, j, k;
    uint32_t numEighVal;
    MathPoint3D point;
    MathVector3D vector;
    std::vector<MathVector3D> vectors;
    std::vector<double> lenghts;
    ListMicroTetraedr::iterator itBegin, itEnd, it;
    double phase = phaseEdit->text().toDouble();
    double maxLenVector = 0.0, minLenVector = DBL_MAX, len, medLen;
    double delta, blueColor, redColor;
    MathPoint3D begin;
    double dx, dy, dz, lenMin, step;

    if(!engine->isCalculate())
    {
        QMessageBox::warning(this, QString("Ошибка"), QString("Сначала нужно нужно выполнить расчет!"), QMessageBox::Ok);
        return;
    }
    numEighVal = numEighValueEdit->text().toInt();
    if(!engine->calculateBasisKoef(numEighVal))
    {
        QMessageBox::warning(this, QString("Ошибка"), QString("Неверный номер собственного значения!"), QMessageBox::Ok);
        return;
    }
    lenWaveEdit->setText(QString("%1").arg(engine->getLenWave()));
    freqEdit->setText(QString("%1").arg(engine->getFrequency()/1000.0));
    setEnabled(false);
    grapher->clearVectors();
    engine->getIteratorTetraedrs(itBegin, itEnd);
    begin = engine->getBeginModelPoint();
    dx = engine->getEndModelPoint().getX() - begin.getX();
    dy = engine->getEndModelPoint().getY() - begin.getY();
    dz = engine->getEndModelPoint().getZ() - begin.getZ();
    lenMin = dx;
    if(dy < lenMin) lenMin = dy;
    if(dz < lenMin) lenMin = dz;
    step = stepDrawEdit->text().toDouble();
    if(step == 0.0)
    {
        step = lenMin/10.0;
        stepDrawEdit->setText(QString("%1").arg(step));
    }
    medLen = 0.0;
    for(i=0; i<dx/step; i++)
        for(j=0; j<dy/step; j++)
            for(k=0; k<dz/step; k++)
            {
                point = engine->getBeginModelPoint() + MathPoint3D(i*step, j*step, k*step);
                for(it = itBegin; it != itEnd; it++)
                    if(it->isInTetraedr(point)) break;
                if(it != itEnd)
                {
                    vector = engine->getField(point, phase, &len);
                    if(vector.getX() != 0.0 || vector.getY() != 0.0 || vector.getZ() != 0.0)
                    {
                        vectors.push_back(vector);
                        lenghts.push_back(len);
                        medLen = medLen + len;
                    }
                }
                qApp->processEvents();
            }

    medLen = medLen/((double)vectors.size());
    for(i=0; i<vectors.size(); i++)
    {
        len = lenghts[i];
        if(!(len > medLen*10.0))
        {
            if(len > maxLenVector) maxLenVector = len;
            if(len < minLenVector) minLenVector = len;
        }
    }
    delta = maxLenVector - minLenVector;
    for(i=0; i<vectors.size(); i++)
    {
        len = vectors[i].lenght();
        if(len > medLen*10.0) len = maxLenVector;
        redColor = len/delta - minLenVector/delta;
        if(redColor > 1.0) redColor = 1.0;
        blueColor = -len/delta + maxLenVector/delta;
        if(blueColor > 1.0) blueColor = 1.0;
        grapher->addVector(vectors[i].unitVector(), grapher->lenEightVector()*len/maxLenVector, QColor(redColor*255.0, 0.0, blueColor*255.0));
        qApp->processEvents();
    }
    addTextConsol(ACTION_END_DRAW_VECTORS);
    setEnabled(true);
    grapher->setEnable(DRAW);
}

void MainDialog::setModelSlot()
{
    grapher->setEnable(POINTS_TYPE);
    grapher->setEnable(LINES_TYPE);
    if(visibleTetraedr) grapher->setEnable(TETRAEDRS_TYPE);
    else grapher->setDisable(TETRAEDRS_TYPE);
    grapher->setDisable(VECTOR_TYPES);
    grapher->setEnable(DRAW);
}

void MainDialog::setVectorSlot()
{
    grapher->setDisable(POINTS_TYPE);
    grapher->setDisable(LINES_TYPE);
    grapher->setDisable(TETRAEDRS_TYPE);
    grapher->setEnable(VECTOR_TYPES);
    grapher->setEnable(DRAW);
}

void MainDialog::setAllSlot()
{
    grapher->setEnable(POINTS_TYPE);
    grapher->setEnable(LINES_TYPE);
    if(visibleTetraedr) grapher->setEnable(TETRAEDRS_TYPE);
    else grapher->setDisable(TETRAEDRS_TYPE);
    grapher->setEnable(VECTOR_TYPES);
    grapher->setEnable(DRAW);
}

void MainDialog::switchTetraedrSlot()
{
    if(visibleTetraedr)
    {
        grapher->setDisable(TETRAEDRS_TYPE);
        grapher->setEnable(DRAW);
        visibleTetraedr = false;
        switchTetraedrButton->setText(QString("Тетраэдры видимы: выкл"));
    }
    else
    {
        grapher->setEnable(TETRAEDRS_TYPE);
        grapher->setEnable(DRAW);
        visibleTetraedr = true;
        switchTetraedrButton->setText(QString("Тетраэдры видимы: вкл"));
    }
}

void MainDialog::saveImageSlot()
{
    QString fileName = QFileDialog::getSaveFileName(this, QString("Сохранить изображение..."));

    if(!grapher->saveImage(fileName))
        QMessageBox::warning(this, QString("Ошибка!"), QString("Файл %1 не сохранен").arg(fileName), QMessageBox::Ok);
}

void MainDialog::addTextConsol(uint32_t numAction)
{
    switch(numAction)
    {
    case ACTION_START_GENERATE:
        consolOut->clear();
        consolOut->insertPlainText(QString("Старт генерации сетки...\n"));
        break;

    case ACTION_START_OBJECT_GENERATE:
        consolOut->insertPlainText(QString("Генерация объектов...\n"));
        break;

    case ACTION_START_PORT_GENERATE:
        consolOut->insertPlainText(QString("Генерация портов...\n"));
        break;

    case ACTION_START_DRAW_OBJECTS:
        consolOut->insertPlainText(QString("Отрисовка объектов в OpenGL...\n"));
        break;

    case ACTION_END_GENERATE:
        consolOut->insertPlainText(QString("Генерация сетки завершена. Для расчета поля нажмите <Рассчитать поле>\n"));
        break;

    case ACTION_START_CALCULATE:
        consolOut->insertPlainText(QString("Старт расчета поля...\n"));
        break;

    case ACTION_START_CREATE_MATRIX:
        consolOut->insertPlainText(QString("Формирование элементов матрицы...\n"));
        break;

    case ACTION_START_SOLVE_MATRIX:
        consolOut->insertPlainText(QString("Решение СЛАУ...\n"));
        break;

    case ACTION_END_CALCULATE:
        consolOut->insertPlainText(QString("Расчет поля завершен. Для отрисовки векторов нажмите <Отрисовать векторы>\n"));
        break;

    case ACTION_END_DRAW_VECTORS:
        consolOut->insertPlainText(QString("Отрисовка векторов завершена.\n"));
        break;

    case ACTION_END_CALCULATE_ERROR:
        consolOut->insertPlainText(QString("Ошибка! Не удалось расчитать поле.\n"));
    }

    consolOut->moveCursor(QTextCursor::End);
}

//EOF
