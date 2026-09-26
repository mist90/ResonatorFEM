#include "MicroEngine.h"
#include <stdio.h>
#include <cstdint>
#include <cmath>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/Dense>

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
    pointGenerate = false;
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
    /* Every boundary node becomes PEC. */
    grid.getIteratorNodes(itBegin, itEnd);
    for(it = itBegin; it != itEnd; ++it)
        if(it->isFlags(NODE_IS_BOUNDARY))
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

void MicroEngine::setResonatorMode(bool enable)
{
    clear();
    resonatorModeEnabled = enable;
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

        MathVector3D c((it->nodes(0)->point().getX() + it->nodes(1)->point().getX() +
                       it->nodes(2)->point().getX() + it->nodes(3)->point().getX()) / 4.0,
                      (it->nodes(0)->point().getY() + it->nodes(1)->point().getY() +
                       it->nodes(2)->point().getY() + it->nodes(3)->point().getY()) / 4.0,
                      (it->nodes(0)->point().getZ() + it->nodes(1)->point().getZ() +
                       it->nodes(2)->point().getZ() + it->nodes(3)->point().getZ()) / 4.0);
        MathVector3D f(0.0, 0.0, 0.0);
        for(i = 0; i < 6; i++)
            f = f + it->getValueBasisFunc(c, i) * basisKoef[6 * it->getSerialNumber() + i];
        points.push_back({ c.getX(), c.getY(), c.getZ() });
        vectors.push_back({ f.getX(), f.getY(), f.getZ() });
    }
    return true;
}

void MicroEngine::clear()
{
    grid.clear();
    globalTable.clear();
    rootsGlobalMatrix.clear();
    basisKoef.clear();
    tableMetallEdges.clear();
    eighValue.clear();
    widthGlobalMatrix = 0;
    pointGenerate = false;
}

bool MicroEngine::_calculate()
{
    /* All tetrahedra are vacuum (relative permittivity 1); the mesh carries no
     * material regions in the current pipeline. */
    std::vector<double> epsilonValuesTetraedrs;

    if(!pointGenerate) return false;
    basisKoef.clear();
    emit startCreateMatrix();
    epsilonValuesTetraedrs.resize(grid.getNumTetraedrs(), 1.0);
    rootsGlobalMatrix.clear();
    globalTable.clear();
    eighValue.clear();
    widthGlobalMatrix = 0;
    calculateMetallEdges();
    return calculateResonatorMode(epsilonValuesTetraedrs);
}

bool MicroEngine::calculateResonatorMode(std::vector<double> &epsilonValuesTetraedrs)
{
    ListMicroTetraedr::iterator itTet, itTetBegin, itTetEnd;
    uint32_t i, j;
    /* local element matrices (6 edges x 6 edges) */
    Eigen::Matrix<double, 6, 6> localMatrixT, localMatrixR;
    std::set<uint32_t> indexes;
    MicroEdge localEdge;
    std::vector<uint32_t> validIndexes;
    uint32_t numEdge;
    uint32_t beginIndexTable;
    uint32_t numLine, numElement;
    uint32_t indexI, indexJ;
    /* для работы с глобальной матрицей */
    uint32_t globalIndex;
    double sign, signLine;
    std::vector<Eigen::Triplet<double> > tripT, tripR;

    /* Assemble the global stiffness (T) and mass (R) directly as Eigen sparse
     * triplets; setFromTriplets sums duplicate (row,col) contributions. */
    grid.getIteratorTetraedrs(itTetBegin, itTetEnd);
    globalIndex = 0;
    uint32_t tetCount = 0;
    for(itTet = itTetBegin; itTet != itTetEnd; itTet++)     /* проход по всем тетраэдрам */
    {
        if((++tetCount & 0x3FF) == 0 && aborted()) return false;   /* cooperative abort */
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
                tripT.emplace_back((int)numLine, (int)numElement,
                                   localMatrixT(validIndexes[i], validIndexes[j]) * sign * signLine);
                tripR.emplace_back((int)numLine, (int)numElement,
                                   localMatrixR(validIndexes[i], validIndexes[j]) * sign * signLine);
            }
        }

    }
    widthGlobalMatrix = globalIndex;
    Eigen::SparseMatrix<double> globalMatrixT((int)globalIndex, (int)globalIndex);
    Eigen::SparseMatrix<double> globalMatrixR((int)globalIndex, (int)globalIndex);
    globalMatrixT.setFromTriplets(tripT.begin(), tripT.end());
    globalMatrixR.setFromTriplets(tripR.begin(), tripR.end());

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
        /* Penalty magnitude must scale like an eigenvalue k^2 = xᵀSx/xᵀRx, i.e.
         * like S/R — NOT like S alone. meanDiag(S) by itself scales as the
         * stiffness (∝ length with the length-scaled Whitney basis), while k^2
         * scales as 1/length^2, so a factor tuned at metre scale collapsed by
         * ~1/length^3 on millimetre geometries and stopped lifting the gradient
         * null space (physical modes then hid behind un-lifted spurious modes).
         * Normalising by meanDiag(R) makes penaltyS ~ mean(k^2), scale-invariant. */
        double meanDiagS = globalMatrixT.diagonal().sum() / (double)widthGlobalMatrix;
        double meanDiagR = globalMatrixR.diagonal().sum() / (double)widthGlobalMatrix;
        double k2Scale   = (meanDiagR > 0.0) ? meanDiagS / meanDiagR : meanDiagS;
        penaltyS = penaltyFactor * k2Scale;    /* lifts null-space modes to ~penaltyS (in k^2 units) */
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
    basisKoef.resize(globalTable.size(), 0.0);
    for(i=0; i<globalTable.size(); i++)
    {
        if(globalTable[i].isFlags(EDGE_NULL)) basisKoef[i] = 0;
        else
        {
            if(globalTable[i].isFlags(EDGE_REVERSE)) basisKoef[i] = -rootsGlobalMatrix[globalTable[i].getNumElement() + numEighVal*widthGlobalMatrix];
            else basisKoef[i] = rootsGlobalMatrix[globalTable[i].getNumElement() + numEighVal*widthGlobalMatrix];
        }
    }
    return true;
}

bool MicroEngine::getEighValues(std::vector<double> &eighValues)
{
    if(!resonatorModeEnabled) return false;
    if(eighValue.size() == 0) return false;
    eighValues = eighValue;
    return true;
}

void MicroEngine::getLocalMatrixResonator(Eigen::Matrix<double, 6, 6> &localMatrixT,
                                          Eigen::Matrix<double, 6, 6> &localMatrixR,
                                          MicroTetraedr &tetraedr,
                                          double epsilon,
                                          std::set<uint32_t> &indexNullElements)
{
    uint32_t i, j;
    /* barycentric-function coefficients */
    Eigen::Matrix4d matrix;
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
        matrix(0, i) = tetraedr.nodes(i)->point().getX();
        matrix(1, i) = tetraedr.nodes(i)->point().getY();
        matrix(2, i) = tetraedr.nodes(i)->point().getZ();
        matrix(3, i) = 1.0;
    }
    volumeTetraedr = fabs(matrix.determinant())/6.0;
    Eigen::Matrix4d minv = matrix.inverse();
    for(i=0; i<4; i++)
    {
        b[i] = minv(i, 0);
        c[i] = minv(i, 1);
        d[i] = minv(i, 2);
    }
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
    localMatrixT.setZero();
    localMatrixR.setZero();
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
            localMatrixT(i, j) = integralD;
            localMatrixR(i, j) = epsilon*integralG;
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

