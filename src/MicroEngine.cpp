#include "MicroEngine.h"
#include <stdio.h>
#include <cstdint>
#include <cmath>
#include <vector>

/* вычисление C для абсорбционных условий */
double getValueC(double *b, double *c, double *d, uint32_t indexPi, uint32_t indexQj, const MathVector3D& normVector)
{
    double cosAlpha = ((MathVector3D)normVector)*MathVector3D(1.0, 0.0, 0.0);
    double cosBeta = ((MathVector3D)normVector)*MathVector3D(0.0, 1.0, 0.0);
    double cosGamma = ((MathVector3D)normVector)*MathVector3D(0.0, 0.0, 1.0);

    return (d[indexPi]*cosBeta - c[indexPi]*cosGamma)*(d[indexQj]*cosBeta - c[indexQj]*cosGamma) +
           (b[indexPi]*cosGamma - d[indexPi]*cosAlpha)*(b[indexQj]*cosGamma - d[indexQj]*cosAlpha) +
           (c[indexPi]*cosAlpha - b[indexPi]*cosBeta)*(c[indexQj]*cosAlpha - b[indexQj]*cosBeta);
}

bool findEdge(std::vector<GlobalTableElement>& table, MicroEdge& edge, uint32_t& getNumEdge)
{
    uint32_t i;

    for(i=0; i<table.size(); i++)
        if(table[i].getEdge() == edge)
        {
            getNumEdge = i;
            return true;
        }
    return false;
}

GlobalTableElement::GlobalTableElement(const MicroEdge &edge, uint32_t numElement, uint32_t flags)
{
    _edge = edge;
    _numElement = numElement;
    _flags = flags;
}

MicroEdge GlobalTableElement::getEdge()
{
    return _edge;
}

uint32_t GlobalTableElement::getNumElement()
{
    return _numElement;
}

void GlobalTableElement::addFlags(uint32_t flags)
{
    _flags |= flags;
}

bool GlobalTableElement::isFlags(uint32_t flags)
{
    if((_flags&flags) == flags) return true;
    else return false;
}

MicroEngine::MicroEngine()
{
    grid.setBeginPoint(MathPoint3D(0.0, 0.0, 0.0));
    grid.setSizeGrid(10.0, 10.0, 10.0);
    pointGenerate = false;
    waveNumber = 1.0;
    indexMainSurface = 0;
    numAction = 0;
    resonatorModeEnabled = false;
    widthGlobalMatrix = 0;
    solveSigmaK2 = -1.0;
    solveNev = 20;
    /* Tree-cotree gauging is OFF by default: zeroing tree DOFs is valid for the
     * driven (source) problem but corrupts the generalized eigenproblem (it
     * restricts the trial space to e_tree=0, which the physical modes do not
     * inhabit). Correct null-space removal needs a mass-metric projection /
     * grad-div penalty (see docs/RESULTS.md). Kept behind a flag for study. */
    useGauge = false;
    penaltyFactor = 0.0;
}

bool MicroEngine::generateFromMesh(const FemMesh &mesh)
{
    std::list<MicroNode>::iterator it, itBegin, itEnd;
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    uint32_t i;

    basisKoef.clear();
    if(!grid.loadTetraMesh(mesh)) return false;
    /* Boundary nodes become PEC (no ports in the mesh-driven resonator path). */
    grid.getIteratorNodes(itBegin, itEnd);
    for(it = itBegin; it != itEnd; ++it)
        if(it->isFlags(NODE_IS_BOUNDARY) && !isPointInEndPort(it->point()))
            it->addFlags(NODE_IS_METALL);
    /* Boundary surfaces become PEC walls. */
    grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
    for(itTet = itTetBegin; itTet != itTetEnd; ++itTet)
        for(i = 0; i < 4; i++)
            if(itTet->isBoundarySurface(i))
                itTet->addFlags(i, SURFACE_IS_METALL);
    calculateMetallEdges();
    if(useGauge) buildTreeCotreeGauge();
    pointGenerate = true;
    return true;
}

bool MicroEngine::calculateSync()
{
    return _calculate();
}

void MicroEngine::setResonatorSolveTarget(double sigmaK2, int nev)
{
    solveSigmaK2 = sigmaK2;
    solveNev = nev;
}

void MicroEngine::setGauge(bool enable)
{
    useGauge = enable;
}

void MicroEngine::setGradDivPenalty(double factor)
{
    penaltyFactor = factor;
}

/* Union-find helpers for spanning-tree construction. */
static uint32_t ufFind(std::vector<uint32_t>& p, uint32_t x)
{
    while(p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}
static void ufUnite(std::vector<uint32_t>& p, uint32_t a, uint32_t b)
{
    uint32_t ra = ufFind(p, a), rb = ufFind(p, b);
    if(ra != rb) p[ra] = rb;
}
static std::pair<uint32_t, uint32_t> edgeSerialKey(MicroEdge e)
{
    uint32_t a = e.getNode1()->getSerialNumber();
    uint32_t b = e.getNode2()->getSerialNumber();
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

void MicroEngine::buildTreeCotreeGauge()
{
    gaugeEdges.clear();
    uint32_t N = grid.getNumNodes();
    if(N == 0) return;

    std::vector<uint32_t> parent(N);
    for(uint32_t i = 0; i < N; i++) parent[i] = i;

    /* PEC edges, keyed by node serials. */
    std::set<std::pair<uint32_t, uint32_t> > metalSet;
    for(uint32_t k = 0; k < tableMetallEdges.size(); k++)
        metalSet.insert(edgeSerialKey(tableMetallEdges[k]));

    /* All unique mesh edges. */
    std::set<std::pair<uint32_t, uint32_t> > allEdges;
    ListMicroTetraedr::iterator it, itBegin, itEnd;
    grid.getIteratorTetraedrs(itBegin, itEnd);
    for(it = itBegin; it != itEnd; ++it)
        for(uint32_t i = 0; i < 6; i++)
            allEdges.insert(edgeSerialKey(it->getEdge(i)));

    /* Pass 1: union all PEC edges first (grounds each conductor; distinct
     * conductors stay in distinct components, preserving inter-conductor modes). */
    std::set<std::pair<uint32_t, uint32_t> >::iterator e;
    for(e = allEdges.begin(); e != allEdges.end(); ++e)
        if(metalSet.find(*e) != metalSet.end())
            ufUnite(parent, e->first, e->second);

    /* Pass 2: spanning tree over non-PEC edges. A non-PEC edge that connects two
     * components is a tree edge -> gauge it to zero; otherwise it is a cotree
     * edge -> keep as a free DOF (carries the physical, circulating field). */
    for(e = allEdges.begin(); e != allEdges.end(); ++e)
        if(metalSet.find(*e) == metalSet.end())
        {
            if(ufFind(parent, e->first) != ufFind(parent, e->second))
            {
                ufUnite(parent, e->first, e->second);
                gaugeEdges.insert(*e);
            }
        }
    fprintf(stderr, "gauge: nodes=%u edges=%zu metal=%zu gauged(tree)=%zu cotree(free)=%zu\n",
            N, allEdges.size(), metalSet.size(), gaugeEdges.size(),
            allEdges.size() - metalSet.size() - gaugeEdges.size());
}

bool MicroEngine::isGaugeEdge(const MicroEdge& edge)
{
    if(!useGauge || gaugeEdges.empty()) return false;
    MicroEdge e = edge;
    return gaugeEdges.find(edgeSerialKey(e)) != gaugeEdges.end();
}

MicroEngine::~MicroEngine()
{
    clear();
}

void MicroEngine::setBeginPoint(const MathPoint3D &point)
{
    grid.setBeginPoint(point);
}

MathPoint3D MicroEngine::getBeginPoint()
{
    return grid.getBeginPoint();
}

MathPoint3D MicroEngine::getBeginModelPoint()
{
    return grid.getBeginModelPoint();
}

MathPoint3D MicroEngine::getEndModelPoint()
{
    return grid.getEndModelPoint();
}

bool MicroEngine::getIteratorTetraedrs(ListMicroTetraedr::iterator &itBegin, ListMicroTetraedr::iterator &itEnd)
{
    return grid.getIteratorTetraedrs(itBegin, itEnd);
}

void MicroEngine::setLenWave(double Value)
{
    if(!resonatorModeEnabled)
        waveNumber = 2.0*M_PI/Value;
}

double MicroEngine::getLenWave()
{
    return 2.0*M_PI/waveNumber;
}

void MicroEngine::setFrequency(double MHzValue)
{
    if(!resonatorModeEnabled)
        waveNumber = 2.0*M_PI*(MHzValue*1E+06)/3E+08;
}

double MicroEngine::getFrequency()
{
    return waveNumber/(2.0*M_PI*1E+06)*3E+08;
}

void MicroEngine::setResonatorMode(bool enable)
{
    clear();
    resonatorModeEnabled = enable;
}

void MicroEngine::setSizeGrid(const double &dx, const double &dy, const double &dz)
{
    clear();
    grid.setSizeGrid(dx, dy, dz);
}

double MicroEngine::getDx()
{
    return grid.getSizeX();
}

double MicroEngine::getDy()
{
    return grid.getSizeY();
}

double MicroEngine::getDz()
{
    return grid.getSizeZ();
}

void MicroEngine::addObject(const MathObject &object, bool solid, const double &sigma, const double &epsilon)
{
    objects.push_back(object);
    sigmaValues.push_back(sigma);
    epsilonValues.push_back(epsilon);
    solidObjects.push_back(solid);
    pointGenerate = false;
}

bool MicroEngine::setMainSurface(uint32_t index)
{
    if(index >= objects.size()) return false;
    indexMainSurface = index;
    return true;
}

void MicroEngine::addPort(const MicroPort &port)
{
    if(!resonatorModeEnabled)
        ports.push_back(port);
}

bool MicroEngine::genPoints()
{
    if(isRunning()) return false;
    numAction = ACTION_GENERATE;
    start();
    return true;
}

bool MicroEngine::calculate()
{
    if(isRunning()) return false;
    numAction = ACTION_CALCULATE;
    start();
    return true;
}

bool MicroEngine::isGenerate()
{
    return pointGenerate;
}

bool MicroEngine::isCalculate()
{
    if(basisKoef.size()) return true;
    else return false;
}

MathVector3D MicroEngine::getField(const MathPoint3D &point, const double &phase, double *outAmpl)
{
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    MathVector3D vectorField, maxVectorField;
    MathComplex<double> koef;
    double ampl, beginPhase;
    uint32_t i;

    vectorField = MathVector3D(0.0, 0.0, 0.0, point);
    maxVectorField = MathVector3D(0.0, 0.0, 0.0, point);
    if(!basisKoef.size()) return vectorField;
    grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
    for(itTet = itTetBegin; itTet != itTetEnd; itTet++)
        if(itTet->isInTetraedr(point)) break;
    if(itTet == itTetEnd) return vectorField;
    for(i=0; i<6; i++)
    {
        koef = basisKoef[6*itTet->getSerialNumber() + i];
        if(resonatorModeEnabled) maxVectorField = maxVectorField + itTet->getValueBasisFunc(point, i)*koef.re();
        else
        {
            // TODO: не реализовано
        }
        ampl = koef.absc();
        beginPhase = koef.phase();
        beginPhase += phase;
        koef.setExponentialArg(ampl, beginPhase);
        vectorField = vectorField + itTet->getValueBasisFunc(point, i)*koef.re();
    }
    if(outAmpl) *outAmpl = maxVectorField.lenght();
    return vectorField;
}

bool MicroEngine::sampleFieldAtCentroids(std::vector<std::array<double, 3> > &points,
                                         std::vector<std::array<double, 3> > &vectors)
{
    ListMicroTetraedr::iterator it, itBegin, itEnd;
    uint32_t i;

    points.clear();
    vectors.clear();
    if(!basisKoef.size()) return false;
    grid.getIteratorTetraedrs(itBegin, itEnd);
    for(it = itBegin; it != itEnd; ++it)
    {
        /* Skip low-quality (sliver) tetrahedra: the Whitney-basis field there
         * scales like 1/volume and produces non-physical spikes. Shape quality
         * q = V / Lrms^3 normalised so a regular tet = 1; drop q < 0.1. */
        double vol = fabs(it->volume());
        double sumL2 = 0.0;
        for(i = 0; i < 6; i++) { double l = it->getEdge(i).lenEdge(); sumL2 += l * l; }
        double lrms3 = pow(sumL2 / 6.0, 1.5);
        double quality = (lrms3 > 0.0) ? (vol / lrms3 / 0.117851) : 0.0;
        if(quality < 0.1) continue;

        MathPoint3D c((it->nodes(0)->point().getX() + it->nodes(1)->point().getX() +
                       it->nodes(2)->point().getX() + it->nodes(3)->point().getX()) / 4.0,
                      (it->nodes(0)->point().getY() + it->nodes(1)->point().getY() +
                       it->nodes(2)->point().getY() + it->nodes(3)->point().getY()) / 4.0,
                      (it->nodes(0)->point().getZ() + it->nodes(1)->point().getZ() +
                       it->nodes(2)->point().getZ() + it->nodes(3)->point().getZ()) / 4.0);
        MathVector3D f(0.0, 0.0, 0.0, c);
        for(i = 0; i < 6; i++)
            f = f + it->getValueBasisFunc(c, i) * basisKoef[6 * it->getSerialNumber() + i].re();
        points.push_back({ c.getX(), c.getY(), c.getZ() });
        vectors.push_back({ f.getX(), f.getY(), f.getZ() });
    }
    return true;
}

void MicroEngine::clear()
{
    grid.clear();
    objects.clear();
    sigmaValues.clear();
    epsilonValues.clear();
    solidObjects.clear();
    ports.clear();
    globalTable.clear();
    rootsGlobalMatrix.clear();
    basisKoef.clear();
    tableMetallEdges.clear();
    eighValue.clear();
    widthGlobalMatrix = 0;
    pointGenerate = false;
    indexMainSurface = 0;
    numAction = 0;
}

void MicroEngine::run()
{
    switch(numAction)
    {
    case ACTION_GENERATE:
        _genPoints();
        emit endGenerate();
        break;

    case ACTION_CALCULATE:
        if(_calculate()) emit endCalculate(true);
        else emit endCalculate(false);
        break;
    }
}

bool MicroEngine::isInObjects(const MathPoint3D &point, const uint32_t indexBegin, bool ignoreMetallAtribut)
{
    uint32_t i;

    if(ignoreMetallAtribut)
    {
        for(i=indexBegin; i<objects.size(); i++)
            if(objects[i].isPointInFigure(point)) return true;
    }
    else
    {
        for(i=indexBegin; i<objects.size(); i++)
            if(objects[i].isPointInFigure(point) && sigmaValues[i] != SIGMA_IDEAL_METALL) return true;
    }
    return false;
}

bool MicroEngine::eraseObject(uint32_t index)
{
    ListMicroTetraedr::iterator itT, itBeginT, itEndT;
    MathPoint3D point;
    bool findOk = false;

    grid.getIteratorTetraedrs(itBeginT, itEndT);
    itT = itBeginT;
    while(itT != itEndT)
    {
        point = (itT->nodes(0)->point() + itT->nodes(1)->point() + itT->nodes(2)->point() + itT->nodes(3)->point())/4.0;
        if(objects[index].isPointInFigure(point) &&
                !isInObjects(point, index + 1, false))
        {
            grid.deleteOneTetraedr(itT);
            findOk = true;
        }
        else itT++;
    }
    return findOk;
}

bool MicroEngine::eraseMetallTetraedrs()
{
    ListMicroTetraedr::iterator itT, itBeginT, itEndT;
    MathPoint3D point;
    bool inPort;
    uint32_t i;

    grid.clearLinkNodes();
    for(i=0; i<objects.size(); i++)
        if(sigmaValues[i] == SIGMA_IDEAL_METALL) eraseObject(i);
    grid.getIteratorTetraedrs(itBeginT, itEndT);
    itT = itBeginT;
    while(itT != itEndT)
    {
        point = (itT->nodes(0)->point() + itT->nodes(1)->point() + itT->nodes(2)->point() + itT->nodes(3)->point())/4.0;
        inPort = false;
        for(i=0; i<ports.size(); i++)
            if(ports[i].isPointInPort(point))
            {
                inPort = true;
                break;
            }
        if(!objects[indexMainSurface].isPointInFigure(point) && !inPort)
        {
            grid.deleteOneTetraedr(itT);
        }
        else itT++;
    }
    return true;
}

void MicroEngine::_genPoints()
{
    uint32_t i, j, k, numBegin;
    std::vector<MathPoint3D> vectorPoints;
    std::list<MicroNode>::iterator it, itBegin, itEnd;
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    MicroNode node;
    bool pointInPort;

    basisKoef.clear();
    grid.clear();
    grid.makeSuperStruct();
    emit startObjectsGenerate();
    /* генерация узлов */
    for(i=0; i<objects.size(); i++)
    {
        vectorPoints.clear();
        node.clear();
        /* добавление нового объекта */
        objects[i].getPoints(vectorPoints, solidObjects[i]);
        for(j=0; j<vectorPoints.size(); j++)
        {
            pointInPort = false;
            if(objects[indexMainSurface].isPointInFigure(vectorPoints[j]))     /* точки не должны выходить за пределы общей поверхности */
            {
                for(k=0; k<ports.size(); k++)                   /* точки не должны попадать в порт */
                    if(ports[k].isPointInPort(vectorPoints[j]))
                    {
                        pointInPort = true;
                        break;
                    }
            }
            else pointInPort = true;
            if(pointInPort) continue;
            if(sigmaValues[i] == SIGMA_IDEAL_METALL) node.addFlags(NODE_IS_METALL);
            node.setPoint(vectorPoints[j]);
            grid.addNode(node);
        }
    }
    emit startPortsGenerate();
    numBegin = objects.size();
    /* генерация объектов портов */
    for(i=0; i<ports.size(); i++)
    {
        ports[i].getObjects(objects, sigmaValues, epsilonValues);
    }
    /* генерация узлов портов */
    for(i=numBegin; i<objects.size(); i++)
    {
        vectorPoints.clear();
        node.clear();

        /* добавление нового объекта */
        objects[i].getPoints(vectorPoints, false);
        for(j=0; j<vectorPoints.size(); j++)
        {
            if(sigmaValues[i] == SIGMA_IDEAL_METALL) node.addFlags(NODE_IS_METALL);
            node.setPoint(vectorPoints[j]);
            grid.addNode(node);
        }
    }

    grid.deleteSuperStruct();       /* удаление суперструктуры */
    eraseMetallTetraedrs();         /* удаление металлических тетраэдров */
    grid.linkNodes();               /* связывание узлов */
    grid.setNumberNodesTetraedrs(); /* последовательная нумерация узлов и тераэдров */
    grid.getIteratorNodes(itBegin, itEnd);
    for(it = itBegin; it != itEnd; it++)
        if(it->isFlags(NODE_IS_BOUNDARY) && !isPointInEndPort(it->point())) it->addFlags(NODE_IS_METALL);
    grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
    for(itTet = itTetBegin; itTet != itTetEnd; itTet++) itTet->calculateBoundariesSurface();
    for(itTet = itTetBegin; itTet != itTetEnd; itTet++)
        for(i=0; i<4; i++)
            if(itTet->isBoundarySurface(i)) itTet->addFlags(i, SURFACE_IS_METALL);
    calculateMetallEdges();
    pointGenerate = true;
}

bool MicroEngine::_calculate()
{
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    uint32_t i;
    /* для получения значений проводимости и диэлектрической проницаемости */
    std::vector<double> sigmaValuesTetraedrs;
    std::vector<double> epsilonValuesTetraedrs;

    if(!pointGenerate) return false;
    basisKoef.clear();
    /* Заполнение структуры значениями проводимости и диэлектрической проницаемости */
    emit startCreateMatrix();
    sigmaValuesTetraedrs.resize(grid.getNumTetraedrs(), 0.0);
    epsilonValuesTetraedrs.resize(grid.getNumTetraedrs(), 1.0);
    for(i=0; i<objects.size(); i++)
    {
        /* проверка на попадание имеющихся тетраэдров в объекты */
        grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
        for(itTet = itTetBegin; itTet != itTetEnd; itTet++)
            if(objects[i].isPointInFigure(itTet->nodes(0)->point()) &&
               objects[i].isPointInFigure(itTet->nodes(1)->point()) &&
               objects[i].isPointInFigure(itTet->nodes(2)->point()) &&
               objects[i].isPointInFigure(itTet->nodes(3)->point()))
            {
                sigmaValuesTetraedrs[itTet->getSerialNumber()] = sigmaValues[i];
                epsilonValuesTetraedrs[itTet->getSerialNumber()] = epsilonValues[i];
            }
    }
    /* очистка массивов */
    rootsGlobalMatrix.clear();
    globalTable.clear();
    eighValue.clear();
    widthGlobalMatrix = 0;
    /* генерация массива металлических тетраэдров */
    calculateMetallEdges();
    /* расчет */
    if(!resonatorModeEnabled) return calculateActiveMode(epsilonValuesTetraedrs, sigmaValuesTetraedrs);
    else return calculateResonatorMode(epsilonValuesTetraedrs);

}

bool MicroEngine::calculateActiveMode(std::vector<double> &epsilonValuesTetraedrs, std::vector<double> &sigmaValuesTetraedrs)
{
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    uint32_t i, j;
    /* для получения локальной матрицы */
    MathMatrix<MathComplex<double> > localMatrix;
    std::set<uint32_t> indexes;
    MicroEdge localEdge;
    std::vector<uint32_t> validIndexes;
    uint32_t numEdge;
    uint32_t beginIndexTable;
    uint32_t numLine, numElement;
    uint32_t indexI, indexJ;
    /* для работы с глобальной матрицей */
    std::vector<GlobalTableElement> globalTable;
    std::vector<MathComplex<double> > rootsGlobalMatrix;
    uint32_t globalIndex;
    MathMatrixSparse<MathComplex<double> > globalMatrix;
    std::vector<MathComplex<double> > globalVector;
    MathComplex<double> valueElement;
    double sign, signLine;

    /* Заполнение глобальной матрицы */
    grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
    globalMatrix.setSize(0, 0);
    globalVector.clear();
    globalIndex = 0;
    for(itTet = itTetBegin; itTet != itTetEnd; itTet++)     /* проход по всем тетраэдрам */
    {
        indexes.clear();
        validIndexes.clear();
        getLocalMatrix(localMatrix, *itTet, sigmaValuesTetraedrs[itTet->getSerialNumber()], epsilonValuesTetraedrs[itTet->getSerialNumber()], indexes);
        /* проход по ребрам тетраэдра для заполнения таблицы соответствия и индексов неметаллических элементов */
        for(i=0; i<6; i++)
        {
            localEdge = itTet->getEdge(i);
            if(indexes.find(i) != indexes.end())    /* если ребро металлическое */
                globalTable.push_back(GlobalTableElement(localEdge, 0, EDGE_NULL));
            else
            {
                validIndexes.push_back(i);
                if(findEdge(globalTable, localEdge, numEdge))
                {
                    if(localEdge.isReverse(globalTable[numEdge].getEdge()))
                        globalTable.push_back(GlobalTableElement(localEdge, globalTable[numEdge].getNumElement(), EDGE_REVERSE));
                    else globalTable.push_back(GlobalTableElement(localEdge, globalTable[numEdge].getNumElement(), 0));
                }
                else
                {
                    globalTable.push_back(GlobalTableElement(localEdge, globalIndex, EDGE_FREE));
                    globalIndex++;
                }
            }
        }
        beginIndexTable = globalTable.size() - 6;
        /* внесение строк в глобальную матрицу */
        for(i=0; i<validIndexes.size(); i++)
        {
            indexI = beginIndexTable + validIndexes[i];
            numLine = globalTable[indexI].getNumElement();
            if(globalTable[indexI].isFlags(EDGE_REVERSE)) signLine = -1.0;
            else signLine = 1.0;
            for(j=0; j<validIndexes.size(); j++)    /* копирование строки из локальной матрицы в глобальную */
            {
                indexJ = beginIndexTable + validIndexes[j];
                numElement = globalTable[indexJ].getNumElement();
                if(globalTable[indexJ].isFlags(EDGE_REVERSE)) sign = -1.0;
                else sign = 1.0;
                if(globalTable[indexI].isFlags(EDGE_FREE) || globalTable[indexJ].isFlags(EDGE_FREE)) valueElement = 0;
                else valueElement = globalMatrix.element(numElement, numLine);
                globalMatrix.setElementExt(numElement, numLine, valueElement + localMatrix.element(validIndexes[j], validIndexes[i])*sign*signLine);
            }
            if(numLine >= globalVector.size()) globalVector.resize(numLine + 1, MathComplex<double>(0.0, 0.0));
            globalVector[numLine] += localMatrix.element(6, validIndexes[i])*signLine;
        }

    }
    /* копирование столбца значений в расширенную матрицу */
    globalMatrix.extensionWidth(1);
    for(i=0; i<globalMatrix.height(); i++)
        globalMatrix.setElement(globalMatrix.width() - 1, i, globalVector[i]);

    /* освобождение памяти */
    globalVector.clear();
    epsilonValuesTetraedrs.clear();
    sigmaValuesTetraedrs.clear();
    /* решение СЛАУ */
    emit startSolveMatrix();
    if(!globalMatrix.solveNoCopyMatrix(rootsGlobalMatrix))
    {
        return false;
    }

    globalMatrix.clear();
    /* заполнение таблицы весовых коэффициентов */
    basisKoef.resize(globalTable.size());
    for(i=0; i<globalTable.size(); i++)
    {
        if(globalTable[i].isFlags(EDGE_NULL)) basisKoef[i] = 0;
        else
        {
            if(globalTable[i].isFlags(EDGE_REVERSE)) basisKoef[i] = -rootsGlobalMatrix[globalTable[i].getNumElement()];
            else basisKoef[i] = rootsGlobalMatrix[globalTable[i].getNumElement()];
        }
    }
    return true;
}

bool MicroEngine::calculateResonatorMode(std::vector<double> &epsilonValuesTetraedrs)
{
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    uint32_t i, j;
    /* для получения локальной матрицы */
    MathMatrix<double> localMatrixT, localMatrixR;
    std::set<uint32_t> indexes;
    MicroEdge localEdge;
    std::vector<uint32_t> validIndexes;
    uint32_t numEdge;
    uint32_t beginIndexTable;
    uint32_t numLine, numElement;
    uint32_t indexI, indexJ;
    /* для работы с глобальной матрицей */
    uint32_t globalIndex;
    MathMatrixSparse<double> globalMatrixT, globalMatrixR;
    double valueElementT, valueElementR;
    double sign, signLine;

    /* Заполнение глобальной матрицы */
    grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
    globalMatrixT.setSize(0, 0);
    globalMatrixR.setSize(0, 0);
    globalIndex = 0;
    for(itTet = itTetBegin; itTet != itTetEnd; itTet++)     /* проход по всем тетраэдрам */
    {
        indexes.clear();
        validIndexes.clear();
        getLocalMatrixResonator(localMatrixT, localMatrixR, *itTet, epsilonValuesTetraedrs[itTet->getSerialNumber()], indexes);
        /* проход по ребрам тетраэдра для заполнения таблицы соответствия и индексов неметаллических элементов */
        for(i=0; i<6; i++)
        {
            localEdge = itTet->getEdge(i);
            if(indexes.find(i) != indexes.end())    /* если ребро металлическое */
                globalTable.push_back(GlobalTableElement(localEdge, 0, EDGE_NULL));
            else
            {
                validIndexes.push_back(i);
                if(findEdge(globalTable, localEdge, numEdge))
                {
                    if(localEdge.isReverse(globalTable[numEdge].getEdge()))
                        globalTable.push_back(GlobalTableElement(localEdge, globalTable[numEdge].getNumElement(), EDGE_REVERSE));
                    else globalTable.push_back(GlobalTableElement(localEdge, globalTable[numEdge].getNumElement(), 0));
                }
                else
                {
                    globalTable.push_back(GlobalTableElement(localEdge, globalIndex, EDGE_FREE));
                    globalIndex++;
                }
            }
        }
        beginIndexTable = globalTable.size() - 6;
        /* внесение строк в глобальную матрицу */
        for(i=0; i<validIndexes.size(); i++)
        {
            indexI = beginIndexTable + validIndexes[i];
            numLine = globalTable[indexI].getNumElement();
            if(globalTable[indexI].isFlags(EDGE_REVERSE)) signLine = -1.0;
            else signLine = 1.0;
            for(j=0; j<validIndexes.size(); j++)    /* копирование строки из локальной матрицы в глобальную */
            {
                indexJ = beginIndexTable + validIndexes[j];
                numElement = globalTable[indexJ].getNumElement();
                if(globalTable[indexJ].isFlags(EDGE_REVERSE)) sign = -1.0;
                else sign = 1.0;
                if(globalTable[indexI].isFlags(EDGE_FREE) || globalTable[indexJ].isFlags(EDGE_FREE))
                {
                    valueElementT = 0;
                    valueElementR = 0;
                }
                else
                {
                    valueElementT = globalMatrixT.element(numElement, numLine);
                    valueElementR = globalMatrixR.element(numElement, numLine);
                }
                globalMatrixT.setElementExt(numElement, numLine, valueElementT + localMatrixT.element(validIndexes[j], validIndexes[i])*sign*signLine);
                globalMatrixR.setElementExt(numElement, numLine, valueElementR + localMatrixR.element(validIndexes[j], validIndexes[i])*sign*signLine);
            }
        }

    }
    widthGlobalMatrix = globalMatrixT.width();

    /* For the grad-div (mass-metric) penalty, extract the discrete gradient G:
     * the reference-orientation node serials of each free DOF, taken from the
     * EDGE_FREE table entries. The solver forms s*(T G) D^-1 (Gᵀ T). */
    std::vector<std::pair<uint32_t, uint32_t> > dofNodes;
    std::vector<double> dofLen;
    std::vector<char> interiorNode;
    double penaltyS = 0.0;
    if(penaltyFactor > 0.0 && widthGlobalMatrix > 0)
    {
        dofNodes.assign(widthGlobalMatrix, std::make_pair(UINT32_MAX, UINT32_MAX));
        dofLen.assign(widthGlobalMatrix, 1.0);
        for(i = 0; i < globalTable.size(); i++)
            if(!globalTable[i].isFlags(EDGE_NULL) && globalTable[i].isFlags(EDGE_FREE))
            {
                MicroEdge e = globalTable[i].getEdge();
                uint32_t dof = globalTable[i].getNumElement();
                dofNodes[dof] = std::make_pair(e.getNode1()->getSerialNumber(), e.getNode2()->getSerialNumber());
                dofLen[dof] = e.lenEdge();
            }
        /* Interior-node mask (metal/PEC nodes are the gauge datum -> excluded). */
        std::list<MicroNode>::iterator itn, itnBegin, itnEnd;
        interiorNode.assign(grid.getNumNodes(), 1);
        grid.getIteratorNodes(itnBegin, itnEnd);
        for(itn = itnBegin; itn != itnEnd; ++itn)
            if(itn->isFlags(NODE_IS_METALL))
                interiorNode[itn->getSerialNumber()] = 0;
        double meanDiag = 0.0;
        for(uint32_t d = 0; d < widthGlobalMatrix; d++) meanDiag += globalMatrixT.element(d, d);
        meanDiag /= (double)widthGlobalMatrix;
        penaltyS = penaltyFactor * meanDiag;   /* lifts null-space modes to ~penaltyS */
    }

    /* освобождение памяти */
    epsilonValuesTetraedrs.clear();
    /* решение СЛАУ и нахождение собственных значений */
    emit startSolveMatrix();
    /* решение задачи TX=k^2 RX. With the grad-div penalty the gradient null space
     * is lifted, so the physical modes are the smallest eigenvalues; otherwise
     * sparse shift-invert around a target, or the dense full-spectrum fallback. */
    bool solved;
    if(penaltyS > 0.0)
        solved = MathEighValVectorShiftInvertGauged(globalMatrixT, globalMatrixR, dofNodes, dofLen,
                     interiorNode, grid.getNumNodes(), penaltyS, (solveSigmaK2 >= 0.0 ? solveSigmaK2 : 0.0),
                     solveNev, eighValue, rootsGlobalMatrix);
    else if(solveSigmaK2 >= 0.0)
        solved = MathEighValVectorShiftInvert(globalMatrixT, globalMatrixR, solveSigmaK2, solveNev, eighValue, rootsGlobalMatrix);
    else
        solved = MathEighValVector(globalMatrixT, globalMatrixR, eighValue, rootsGlobalMatrix);
    if(!solved)
    {
        return false;
    }
    else
    {
        globalMatrixT.clear();
        globalMatrixR.clear();
        /* заполнение таблицы весовых коэффициентов */
        calculateBasisKoef(0);
    }
    return true;
}

bool MicroEngine::calculateBasisKoef(uint32_t numEighVal)
{
    uint32_t i;

    if(!resonatorModeEnabled) return false;
    if(rootsGlobalMatrix.size() == 0) return false;
    if(numEighVal >= eighValue.size()) return false;
    waveNumber = sqrt(eighValue[numEighVal]);
    basisKoef.resize(globalTable.size(), MathComplex<double>(0.0, 0.0));
    for(i=0; i<globalTable.size(); i++)
    {
        if(globalTable[i].isFlags(EDGE_NULL)) basisKoef[i] = 0;
        else
        {
            if(globalTable[i].isFlags(EDGE_REVERSE)) basisKoef[i] = -rootsGlobalMatrix[globalTable[i].getNumElement() + numEighVal*widthGlobalMatrix];
            else basisKoef[i] = rootsGlobalMatrix[globalTable[i].getNumElement() + numEighVal*widthGlobalMatrix];
        }
    }
    waveNumber = sqrt(eighValue[numEighVal]);
    return true;
}

bool MicroEngine::getEighValues(std::vector<double> &eighValues)
{
    if(!resonatorModeEnabled) return false;
    if(eighValue.size() == 0) return false;
    eighValues = eighValue;
    return true;
}

void MicroEngine::getLocalMatrix(MathMatrix<MathComplex<double> > &localMatrix,
                                 MicroTetraedr& tetraedr,
                                 double sigma,
                                 double epsilon,
                                 std::set<uint32_t> &indexNullElements)
{
    /* TODO: в классе MathMatrix сначала идет столбец, потом строка */
    uint32_t i, j, k;
    /* для коэффициенты барицентрических функций */
    MathMatrix<double> matrix(4, 4);
    double b[4], c[4], d[4];
    /* для вычисления вспомогательных величин */
    MathVector3D addV[4][4];
    double addPhi[4][4];
    double matrixC[4][4];
    /* для вычисления элементов локальной матрицы */
    double integralD, integralG;
    double volumeTetraedr;
    uint32_t m1, m2, n1, n2;
    /* для вычисления свободного столбца и абсорбционных условий */
    uint32_t numEdge, numPoints[3];
    MathVector3D vectorsH[3];
    MathVector3D normVector;
    MathPoint3D pointsField[3];
    MathPoint3D centerPoint;
    double phases[3];
    double areaSurface;
    MathComplex<double> sElement, cmplxPhases[3];
    double valueC;

    /* определение ребер, граничащих с металлом */
    indexNullElements.clear();
    for(i=0; i<6; i++)
        if(isMetallEdge(tetraedr.getEdge(i)))
            indexNullElements.insert(i);

    /* вычисление барицентрических коэффициентов */
    for(i=0; i<4; i++)
    {
        matrix.element(i, 0) = tetraedr.nodes(i)->point().getX();
        matrix.element(i, 1) = tetraedr.nodes(i)->point().getY();
        matrix.element(i, 2) = tetraedr.nodes(i)->point().getZ();
        matrix.element(i, 3) = 1.0;
    }
    volumeTetraedr = fabs(matrix.determinant())/6.0;
    matrix = matrix.inverseMatrix();
    for(i=0; i<4; i++)
    {
        b[i] = matrix.element(0, i);
        c[i] = matrix.element(1, i);
        d[i] = matrix.element(2, i);
    }
    matrix.clear();
    /* вычисление дополнительных величин */
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
        {
            addV[i][j] = MathVector3D(c[i]*d[j] - c[j]*d[i], b[j]*d[i] - b[i]*d[j], b[i]*c[j] - b[j]*c[i]);
            addPhi[i][j] = b[i]*b[j] + c[i]*c[j] + d[i]*d[j];
            matrixC[i][j] = 1.0/20.0;
        }
    for(i=0; i<4; i++) matrixC[i][i] = matrixC[i][i]*2.0;
    /* вычисление интегралов и элементов матрицы */
    localMatrix.setSize(7, 6);
    for(i=0; i<6; i++)
    {
        if(indexNullElements.find(i) != indexNullElements.end()) continue;
        for(j=0; j<6; j++)
        {
            if(indexNullElements.find(j) != indexNullElements.end()) continue;
            m1 = tetraedr.getNumBeginEdge(i);
            m2 = tetraedr.getNumEndEdge(i);
            n1 = tetraedr.getNumBeginEdge(j);
            n2 = tetraedr.getNumEndEdge(j);
            integralD = 4.0*volumeTetraedr*(tetraedr.getEdge(i).lenEdge())*(tetraedr.getEdge(j).lenEdge()) * (addV[m1][m2]*addV[n1][n2]);
            integralG = volumeTetraedr*(tetraedr.getEdge(i).lenEdge())*(tetraedr.getEdge(j).lenEdge())*
                    (addPhi[m2][n2]*matrixC[m1][n1] - addPhi[m2][n1]*matrixC[m1][n2] - addPhi[m1][n2]*matrixC[m2][n1] + addPhi[m1][n1]*matrixC[m2][n2]);
            localMatrix.element(j, i) = MathComplex<double>(integralD) +
                                MathComplex<double>(0.0, waveNumber*120.0*M_PI*sigma*integralG) -
                                MathComplex<double>(waveNumber*waveNumber*epsilon*integralG);
        }
    }
    /* вычисление свободного столбца */
    for(i=0; i<4; i++)              /* проход по плоскостям тетраэдра */
    {
        for(j=0; j<3; j++)          /* получение точек плоскости */
            pointsField[j] = tetraedr.nodes(tetraedr.getNumPointFromSurface(i, j))->point();
        if(getFieldPoints(pointsField[0], vectorsH[0], phases[0]) &&
           getFieldPoints(pointsField[1], vectorsH[1], phases[1]) &&
           getFieldPoints(pointsField[2], vectorsH[2], phases[2]))
        {
            centerPoint = (pointsField[0] + pointsField[1] + pointsField[2])/3.0;
            areaSurface = AreaTriangle(pointsField[0], pointsField[1], pointsField[2]);
            for(j=0; j<3; j++) cmplxPhases[j].setExponentialArg(1.0, phases[j]);
            for(j=0; j<3; j++)
            {
                numEdge = tetraedr.getNumEdgeFromSurface(i, j);
                sElement = (cmplxPhases[0] + cmplxPhases[1] + cmplxPhases[2])/3.0;
                sElement = sElement*(((vectorsH[0] + vectorsH[1] + vectorsH[2])*(1.0/3.0))*tetraedr.getValueBasisFunc(centerPoint, numEdge))*
                        areaSurface*waveNumber*120.0*M_PI;
                localMatrix.element(6, numEdge) = sElement;
            }
            break;
        }
    }
    /* добавление в матрицу слагаемых, задаваемых абсорбционными условиями */
    for(i=0; i<4; i++)              // проход по плоскостям тетраэдра
    {
        for(j=0; j<3; j++)          // получение точек плоскости
        {
            numPoints[j] = tetraedr.getNumPointFromSurface(i, j);
            pointsField[j] = tetraedr.nodes(numPoints[j])->point();
        }
        normVector = tetraedr.getNormVector(i);
        if(getAbsorbPoints(pointsField[0]) &&
                   getAbsorbPoints(pointsField[1]) &&
                   getAbsorbPoints(pointsField[2]))
        {
            areaSurface = AreaTriangle(pointsField[0], pointsField[1], pointsField[2]);
            for(j=0; j<3; j++)
            {
                numEdge = tetraedr.getNumEdgeFromSurface(i, j);
                for(k=0; k<6; k++)
                {
                    m1 = tetraedr.getNumBeginEdge(numEdge);
                    m2 = tetraedr.getNumEndEdge(numEdge);
                    n1 = tetraedr.getNumBeginEdge(k);
                    n2 = tetraedr.getNumEndEdge(k);
                    valueC = getValueC(b, c, d, m2, n2, normVector) - getValueC(b, c, d, m2, n1, normVector) - getValueC(b, c, d, m1, n2, normVector) + getValueC(b, c, d, m1, n1, normVector);
                    localMatrix.element(k, numEdge) = localMatrix.element(k, numEdge) + MathComplex<double>(0.0, waveNumber/12.0*sqrt(epsilon)*areaSurface*valueC);
                }
            }
            break;
        }
    }
}

void MicroEngine::getLocalMatrixResonator(MathMatrix<double> &localMatrixT,
                                          MathMatrix<double> &localMatrixR,
                                          MicroTetraedr &tetraedr,
                                          double epsilon,
                                          std::set<uint32_t> &indexNullElements)
{
    /* TODO: в классе MathMatrix сначала идет столбец, потом строка */
    uint32_t i, j;
    /* для коэффициенты барицентрических функций */
    MathMatrix<double> matrix(4, 4);
    double b[4], c[4], d[4];
    /* для вычисления вспомогательных величин */
    MathVector3D addV[4][4];
    double addPhi[4][4];
    double matrixC[4][4];
    /* для вычисления элементов локальной матрицы */
    double integralD, integralG;
    double volumeTetraedr;
    uint32_t m1, m2, n1, n2;

    /* определение ребер, граничащих с металлом, и калибровочных (tree) ребер */
    indexNullElements.clear();
    for(i=0; i<6; i++)
    {
        MicroEdge edgeI = tetraedr.getEdge(i);
        if(isMetallEdge(edgeI) || isGaugeEdge(edgeI))
            indexNullElements.insert(i);
    }

    /* вычисление барицентрических коэффициентов */
    for(i=0; i<4; i++)
    {
        matrix.element(i, 0) = tetraedr.nodes(i)->point().getX();
        matrix.element(i, 1) = tetraedr.nodes(i)->point().getY();
        matrix.element(i, 2) = tetraedr.nodes(i)->point().getZ();
        matrix.element(i, 3) = 1.0;
    }
    volumeTetraedr = fabs(matrix.determinant())/6.0;
    matrix = matrix.inverseMatrix();
    for(i=0; i<4; i++)
    {
        b[i] = matrix.element(0, i);
        c[i] = matrix.element(1, i);
        d[i] = matrix.element(2, i);
    }
    matrix.clear();
    /* вычисление дополнительных величин */
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
        {
            addV[i][j] = MathVector3D(c[i]*d[j] - c[j]*d[i], b[j]*d[i] - b[i]*d[j], b[i]*c[j] - b[j]*c[i]);
            addPhi[i][j] = b[i]*b[j] + c[i]*c[j] + d[i]*d[j];
            matrixC[i][j] = 1.0/20.0;
        }
    for(i=0; i<4; i++) matrixC[i][i] = matrixC[i][i]*2.0;
    /* вычисление интегралов и элементов матрицы */
    localMatrixT.setSize(6, 6);
    localMatrixR.setSize(6, 6);
    for(i=0; i<6; i++)
    {
        if(indexNullElements.find(i) != indexNullElements.end()) continue;
        for(j=0; j<6; j++)
        {
            if(indexNullElements.find(j) != indexNullElements.end()) continue;
            m1 = tetraedr.getNumBeginEdge(i);
            m2 = tetraedr.getNumEndEdge(i);
            n1 = tetraedr.getNumBeginEdge(j);
            n2 = tetraedr.getNumEndEdge(j);
            integralD = 4.0*volumeTetraedr*(tetraedr.getEdge(i).lenEdge())*(tetraedr.getEdge(j).lenEdge()) * (addV[m1][m2]*addV[n1][n2]);
            integralG = volumeTetraedr*(tetraedr.getEdge(i).lenEdge())*(tetraedr.getEdge(j).lenEdge())*
                    (addPhi[m2][n2]*matrixC[m1][n1] - addPhi[m2][n1]*matrixC[m1][n2] - addPhi[m1][n2]*matrixC[m2][n1] + addPhi[m1][n1]*matrixC[m2][n2]);
            integralG = integralG;
            localMatrixT.element(j, i) = integralD;
            localMatrixR.element(j, i) = epsilon*integralG;
        }
    }
}

bool MicroEngine::calculateMetallEdges()
{
    ListMicroTetraedr::iterator it, itBegin, itEnd;
    MicroEdge edge;
    uint32_t i, j;

    tableMetallEdges.clear();
    grid.getIteratorTetraedrs(itBegin, itEnd);
    for(it=itBegin; it!=itEnd; it++)
    {
        for(i=0; i<4; i++)
            if(it->isFlags(i, SURFACE_IS_METALL))
                for(j=0; j<3; j++)
                {
                    edge = MicroEdge(it->nodes(i + j), it->nodes(i + (j + 1)%3));
                    if(!isPointInEndPort(edge.getNode1()->point()) || !isPointInEndPort(edge.getNode2()->point()))
                        if(!isMetallEdge(edge)) tableMetallEdges.push_back(edge);
                }
    }
    return true;
}

bool MicroEngine::isMetallEdge(const MicroEdge &edge)
{
    uint32_t i;

    for(i=0; i<tableMetallEdges.size(); i++)
        if(tableMetallEdges[i] == edge) return true;
    return false;
}

bool MicroEngine::getFieldPoints(const MathPoint3D &point, MathVector3D &amplVectorH, double &phaseVectorH)
{
    uint32_t i;

    for(i=0; i<ports.size(); i++)
        if(ports[i].getFieldPoints(point, amplVectorH, phaseVectorH)) return true;
    return false;
}


bool MicroEngine::getAbsorbPoints(const MathPoint3D &point)
{
    uint32_t i;

    for(i=0; i<ports.size(); i++)
        if(ports[i].getAbsorbPoints(point)) return true;
    return false;
}

bool MicroEngine::isPointInEndPort(const MathPoint3D &point)
{
    uint32_t i;

    for(i=0; i<ports.size(); i++)
        if(ports[i].isPointInEndPort(point)) return true;
    return false;
}

//EOF
