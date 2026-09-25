#ifndef MICROENGINE_H
#define MICROENGINE_H
#include <vector>
#include <array>
#include <set>
#include <utility>
#include <QObject>
#include "MicroGrid.h"
#include "MicroNode.h"
#include "MathEighValues.h"
#include "FemMesh.h"
#include <Eigen/Core>

#define EDGE_REVERSE        (0x1)
#define EDGE_NULL           (0x2)
#define EDGE_FREE           (0x4)

/* Links a tetrahedron's local edge to a global DOF in the assembled matrices. */
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

/* Edge-element (Whitney) FEM solver for cavity resonator eigenmodes. Driven by an
 * external conforming tetrahedral mesh (from Gmsh); graphics-free. */
class MicroEngine : public QObject
{
    Q_OBJECT
public:
    MicroEngine();
    ~MicroEngine();

    /* Load a conforming tet mesh and flag every exterior (PEC) boundary. */
    bool            generateFromMesh(const FemMesh& mesh);
    /* Run the resonator eigen-solve synchronously on the calling thread. */
    bool            calculateSync();
    void            setResonatorMode(bool enable);

    /* Sparse shift-invert target: `nev` modes with k^2 nearest `sigmaK2`
     * (sigmaK2 < 0 falls back to the dense full-spectrum solver). */
    void            setResonatorSolveTarget(double sigmaK2, int nev);
    /* Grad-div penalty (factor * meanDiag(S)): lifts the gradient null space so the
     * physical modes become the smallest eigenvalues. 0 disables it. */
    void            setGradDivPenalty(double factor);
    /* Tree-cotree gauging — experimental, off by default (see MicroEngine.cpp). */
    void            setGauge(bool enable);

    bool            getEighValues(std::vector<double>& eighValues);
    /* Reconstruct the basis coefficients of eigenmode `numEighVal` (call before
     * sampling its field). */
    bool            calculateBasisKoef(uint32_t numEighVal);
    /* Sample the current mode's real field vector at every tetrahedron centroid
     * (one pass, low-quality slivers skipped) — for the field-glyph view. */
    bool            sampleFieldAtCentroids(std::vector<std::array<double, 3> >& points,
                                           std::vector<std::array<double, 3> >& vectors);
    void            clear();

signals:
    void            startCreateMatrix();
    void            startSolveMatrix();

private:
    bool            _calculate();
    bool            calculateResonatorMode(std::vector<double>& epsilonValuesTetraedrs);
    void            getLocalMatrixResonator(Eigen::Matrix<double, 6, 6>& localMatrixT,
                                            Eigen::Matrix<double, 6, 6>& localMatrixR,
                                            MicroTetraedr& tetraedr,
                                            double epsilon,
                                            std::set<uint32_t>& indexNullElements);
    bool            calculateMetallEdges();
    bool            isMetallEdge(const MicroEdge& edge);
    /* Tree-cotree gauge: spanning tree of the mesh-edge graph (PEC edges first). */
    void            buildTreeCotreeGauge();
    bool            isGaugeEdge(const MicroEdge& edge);

    bool            pointGenerate;
    bool            resonatorModeEnabled;
    MicroGrid       grid;
    /* Assembled-system state. */
    std::vector<GlobalTableElement>   globalTable;
    std::vector<double>               rootsGlobalMatrix;
    std::vector<double>               eighValue;
    uint32_t                          widthGlobalMatrix;
    std::vector<MicroEdge>            tableMetallEdges;
    std::vector<double> basisKoef;   /* real edge weights of the current mode */
    /* Eigensolver target: solveSigmaK2 >= 0 => sparse shift-invert, else dense. */
    double          solveSigmaK2;
    int             solveNev;
    /* Tree-cotree gauge state (gauged non-PEC tree edges, keyed by node serials). */
    bool            useGauge;
    std::set<std::pair<uint32_t, uint32_t> > gaugeEdges;
    /* Grad-div penalty weight factor (0 = off). */
    double          penaltyFactor;
};

#endif // MICROENGINE_H
