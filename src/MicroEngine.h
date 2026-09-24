#ifndef MICROENGINE_H
#define MICROENGINE_H
#include <vector>
#include <set>
#include <QThread>
#include <QMessageBox>
#include "MicroGrid.h"
#include "MicroNode.h"
#include "MathGrapher.h"
#include "MathObject.h"
#include "MicroPort.h"
#include "MathMatrix.h"
#include "MathMatrixSparse.h"
#include "MathEighValues.h"
#include "FemMesh.h"

#define EDGE_REVERSE        (0x1)
#define EDGE_NULL           (0x2)
#define EDGE_FREE           (0x4)

#define ACTION_GENERATE     (0x1)
#define ACTION_CALCULATE    (0x2)
#define ACTION_DRAW         (0x3)

/* элемент для таблицы связывания глобальной и локальных матриц */
class GlobalTableElement
{
public:
    GlobalTableElement(const MicroEdge& edge, uint32_t numElement, uint32_t flags);
    MicroEdge   getEdge();
    uint32_t    getNumElement();
    void        addFlags(uint32_t flags);
    bool        isFlags(uint32_t flags);
private:
    MicroEdge   _edge;
    uint32_t    _numElement;
    uint32_t    _flags;
};

class MicroEngine:public QThread
{
    Q_OBJECT
public:
    MicroEngine();                          /* headless (no grapher) */
    MicroEngine(MathGrapher& grapher);
    ~MicroEngine();
    /* Mesh-driven resonator path (replaces genPoints for externally meshed
     * geometry). generateFromMesh loads a conforming tet mesh and flags PEC
     * boundaries; calculateSync runs the eigen-solve on the calling thread. */
    bool            generateFromMesh(const FemMesh& mesh);
    bool            calculateSync();
    /* Select the sparse shift-invert eigensolver (resonator mode): find `nev`
     * modes with k^2 nearest `sigmaK2`. Targets the physical range directly and
     * scales to fine meshes. Pass sigmaK2 < 0 to fall back to the dense solver. */
    void            setResonatorSolveTarget(double sigmaK2, int nev);
    /* Enable/disable tree-cotree gauging (default off — see MicroEngine.cpp). */
    void            setGauge(bool enable);
    /* Grad-div penalty (resonator mode): add factor * meanDiag(S) * G Gᵀ to the
     * stiffness, where G is the discrete gradient. Lifts the gradient null-space
     * modes to higher k² so the physical modes (nearly divergence-free) become
     * the smallest eigenvalues — needed for very-low-frequency modes (spiral).
     * factor = 0 disables it (default). */
    void            setGradDivPenalty(double factor);
    void            setSizeGrid(const double& dx, const double& dy, const double& dz);
    double          getDx();
    double          getDy();
    double          getDz();
    void            setBeginPoint(const MathPoint3D& point);
    MathPoint3D     getBeginPoint();
    MathPoint3D     getBeginModelPoint();
    MathPoint3D     getEndModelPoint();
    bool            getIteratorTetraedrs(ListMicroTetraedr::iterator& itBegin, ListMicroTetraedr::iterator& itEnd);
    void            setLenWave(double Value);
    double          getLenWave();
    void            setFrequency(double MHzValue);      /* задание частоты в МГц */
    double          getFrequency();                     /* получение частоты в МГц */
    void            setResonatorMode(bool enable);
    bool            getEighValues(std::vector<double>& eighValues);
    void            addObject(const MathObject& object, bool solid, const double& sigma, const double& epsilon);
    bool            drawObject(uint32_t index, QColor color);
    bool            drawAllObject();
    bool            setMainSurface(uint32_t index);
    void            addPort(const MicroPort& port);
    bool            genPoints();            /* генерация точек */
    bool            calculate();            /* расчет поля */
    bool            calculateBasisKoef(uint32_t numEighVal); /* расчет базисных коэффициентов (только резонаторная задача) */
    bool            isGenerate();
    bool            isCalculate();
    MathVector3D    getField(const MathPoint3D& point, const double& phase, double* outAmpl = 0);
    void            clear();
signals:
    void            startObjectsGenerate();
    void            startPortsGenerate();
    void            startDrawObjects();
    void            endGenerate();
    void            startCreateMatrix();
    void            startSolveMatrix();
    void            endCalculate(bool ok);
protected:
    void            run();
private:
    bool            isInObjects(const MathPoint3D& point, const uint32_t indexBegin, bool ignoreMetallAtribut = true);
    bool            eraseObject(uint32_t index);
    bool            eraseMetallTetraedrs();
    void            _genPoints();
    bool            _calculate();
    bool            calculateActiveMode(std::vector<double>& epsilonValuesTetraedrs, std::vector<double>& sigmaValuesTetraedrs);
    bool            calculateResonatorMode(std::vector<double>& epsilonValuesTetraedrs);
    /* Получение локальной матрицы элемента */
    void            getLocalMatrix(MathMatrix<MathComplex<double> >& localMatrix,
                                                    MicroTetraedr& tetraedr,
                                                    double sigma,
                                                    double epsilon,
                                                    std::set<uint32_t>& indexNullElements);
    void            getLocalMatrixResonator(MathMatrix<double>& localMatrixT,
                                                    MathMatrix<double>& localMatrixR,
                                                    MicroTetraedr& tetraedr,
                                                    double epsilon,
                                                    std::set<uint32_t>& indexNullElements);
    bool            calculateMetallEdges();
    bool            isMetallEdge(const MicroEdge& edge);
    /* Tree-cotree gauge: build a spanning tree of the mesh-edge graph (PEC edges
     * first) and record the non-PEC tree edges to be gauged to zero. */
    void            buildTreeCotreeGauge();
    bool            isGaugeEdge(const MicroEdge& edge);
    bool            getFieldPoints(const MathPoint3D &point, MathVector3D &amplVectorH, double &phaseVectorH);
    bool            getAbsorbPoints(const MathPoint3D& point);
    bool            isPointInEndPort(const MathPoint3D& point);
    bool            pointGenerate;
    bool            resonatorModeEnabled;
    uint32_t        numAction;
    MicroGrid       grid;
    MathGrapher     *_grapher;
    uint32_t        indexMainSurface;
    /* Для резонансной задачи */
    std::vector<GlobalTableElement> globalTable;
    std::vector<double> rootsGlobalMatrix;
    std::vector<double> eighValue;
    uint32_t        widthGlobalMatrix;
    /* волновое число */
    double          waveNumber;
    std::vector<MathObject> objects;
    /* Значения диэлектрической проницаемости и проводимости */
    std::vector<double> sigmaValues;
    std::vector<double> epsilonValues;
    std::vector<bool> solidObjects;
    /* Граничные значения вектора E */
    std::vector<MicroPort> ports;
    /* Решения СЛАУ - весовые коэффициенты базисных функций */
    std::vector<MicroEdge> tableMetallEdges;
    std::vector<MathComplex<double> > basisKoef;
    /* Resonator eigensolver target: if solveSigmaK2 >= 0 use sparse shift-invert
     * around solveSigmaK2 for solveNev modes, else dense full-spectrum. */
    double          solveSigmaK2;
    int             solveNev;
    /* Tree-cotree gauging state. gaugeEdges holds gauged (tree) non-PEC edges,
     * keyed by sorted node serial numbers. */
    bool            useGauge;
    std::set<std::pair<uint32_t, uint32_t> > gaugeEdges;
    /* Grad-div penalty weight factor (0 = off). */
    double          penaltyFactor;
};



#endif // MICROENGINE_H
