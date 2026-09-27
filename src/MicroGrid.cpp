#include "MicroGrid.h"
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <cstdio>


MicroGrid::MicroGrid()
{
    isClearTetraedrs = false;
    isMakeSuperStruct = false;
    isLinkNodes = true;
    isNumNodesTetraedrs = false;
}


bool MicroGrid::setNumberNodesTetraedrs()
{
    std::list<MicroNode>::iterator itNode;
    ListMicroTetraedr::iterator itTetraedr;
    uint32_t count = 0;

    if(isNumNodesTetraedrs) return true;
    if(isClearTetraedrs) return false;
    for(itNode = listNodes.begin(); itNode != listNodes.end(); itNode++)
    {
        itNode->_serialNumber = count;
        count++;
    }
    count = 0;
    for(itTetraedr = listTetraedr.begin(); itTetraedr != listTetraedr.end(); itTetraedr++)
    {
        itTetraedr->setSerialNumber(count);
        count++;
    }
    isNumNodesTetraedrs = true;
    return true;
}

bool MicroGrid::getIteratorNodes(std::list<MicroNode>::iterator &itBegin, std::list<MicroNode>::iterator& itEnd)
{
    if(!listNodes.size()) return false;
    itBegin = listNodes.begin();
    itEnd = listNodes.end();
    return true;
}

bool MicroGrid::getIteratorTetraedrs(ListMicroTetraedr::iterator& itBegin, ListMicroTetraedr::iterator& itEnd)
{
    if(!listTetraedr.size()) return false;
    itBegin = listTetraedr.begin();
    itEnd = listTetraedr.end();
    return true;
}

uint32_t MicroGrid::getNumNodes()
{
    return listNodes.size();
}

uint32_t MicroGrid::getNumTetraedrs()
{
    return listTetraedr.size();
}

void MicroGrid::clear()
{
    isClearTetraedrs = false;
    isMakeSuperStruct = false;
    isLinkNodes = false;
    isNumNodesTetraedrs = false;
    listNodes.clear();
    listTetraedr.clear();
}

bool MicroGrid::loadTetraMesh(const FemMesh &mesh)
{
    if(mesh.tets.empty()) return false;
    clear();

    /* 1. Nodes. std::list iterators are stable, so we cache one per node index. */
    std::vector<std::list<MicroNode>::iterator> nodeIt;
    nodeIt.reserve(mesh.nodes.size());
    for(uint32_t i = 0; i < mesh.nodes.size(); i++)
    {
        MicroNode node;
        node.setPoint(MathVector3D(mesh.nodes[i][0], mesh.nodes[i][1], mesh.nodes[i][2]));
        listNodes.push_back(node);
        std::list<MicroNode>::iterator it = listNodes.end();
        --it;
        nodeIt.push_back(it);
    }

    /* 2. Tetrahedra, referencing the cached node iterators. Degenerate (sliver)
     * tetrahedra are skipped: their edges would carry ~0 stiffness and mass,
     * producing zero rows that make the FEM pencil structurally singular (breaks
     * both the dense Cholesky and the sparse shift-invert factorization). The
     * threshold is relative to the largest element so it is scale-independent. */
    double maxVol = 0.0;
    for(uint32_t i = 0; i < mesh.tets.size(); i++)
    {
        const std::array<uint32_t, 4> &t = mesh.tets[i];
        MathVector3D p0 = nodeIt[t[0]]->point(), p1 = nodeIt[t[1]]->point();
        MathVector3D p2 = nodeIt[t[2]]->point(), p3 = nodeIt[t[3]]->point();
        double v = fabs(VolumeTetraedr(p0, p1, p2, p3));
        if(v > maxVol) maxVol = v;
    }
    const double volTol = 1.0e-9 * maxVol;

    std::vector<ListMicroTetraedr::iterator>   tetIt;
    std::vector<std::array<uint32_t, 4> >      keptTets;
    tetIt.reserve(mesh.tets.size());
    keptTets.reserve(mesh.tets.size());
    uint32_t nDegenerate = 0;
    for(uint32_t i = 0; i < mesh.tets.size(); i++)
    {
        const std::array<uint32_t, 4> &t = mesh.tets[i];
        MathVector3D p0 = nodeIt[t[0]]->point(), p1 = nodeIt[t[1]]->point();
        MathVector3D p2 = nodeIt[t[2]]->point(), p3 = nodeIt[t[3]]->point();
        if(fabs(VolumeTetraedr(p0, p1, p2, p3)) <= volTol) { nDegenerate++; continue; }
        MicroTetraedr tet;
        for(uint32_t k = 0; k < 4; k++)
            tet.nodes(k) = nodeIt[t[k]];
        listTetraedr.push_back(tet);
        ListMicroTetraedr::iterator it = listTetraedr.end();
        --it;
        tetIt.push_back(it);
        keptTets.push_back(t);
    }
    if(nDegenerate) fprintf(stderr, "loadTetraMesh: skipped %u degenerate tets (of %zu)\n",
                            nDegenerate, mesh.tets.size());

    /* 3. Face adjacency. Surface i of a tet is the face on local nodes
     * {i, (i+1)%4, (i+2)%4} (the convention used by calculateBoundariesSurface
     * and MicroEngine::calculateMetallEdges). A face key is the sorted triple of
     * global node indices; a face seen twice links two tetrahedra, a face seen
     * once is a boundary. We erase matched faces so the leftovers are boundaries. */
    std::map<std::array<uint32_t, 3>, std::pair<uint32_t, uint32_t> > faceMap;
    for(uint32_t ti = 0; ti < keptTets.size(); ti++)
    {
        const std::array<uint32_t, 4> &t = keptTets[ti];
        for(uint32_t i = 0; i < 4; i++)
        {
            std::array<uint32_t, 3> key = { t[i], t[(i + 1) % 4], t[(i + 2) % 4] };
            std::sort(key.begin(), key.end());
            std::map<std::array<uint32_t, 3>, std::pair<uint32_t, uint32_t> >::iterator f = faceMap.find(key);
            if(f == faceMap.end())
                faceMap[key] = std::make_pair(ti, i);
            else
            {
                uint32_t tj = f->second.first;
                uint32_t sj = f->second.second;
                tetIt[ti]->tetraedrs(i)  = tetIt[tj];
                tetIt[tj]->tetraedrs(sj) = tetIt[ti];
                faceMap.erase(f);
            }
        }
    }

    /* 4. Remaining entries are boundary faces: flag their nodes as boundary. */
    for(std::map<std::array<uint32_t, 3>, std::pair<uint32_t, uint32_t> >::iterator f = faceMap.begin();
        f != faceMap.end(); ++f)
        for(int k = 0; k < 3; k++)
            nodeIt[f->first[k]]->addFlags(NODE_IS_BOUNDARY);

    /* 5. Sequential node/tetrahedron numbering. */
    setNumberNodesTetraedrs();

    /* 6. Mark boundary surfaces (needs NODE_IS_BOUNDARY set above). */
    for(ListMicroTetraedr::iterator it = listTetraedr.begin(); it != listTetraedr.end(); ++it)
        it->calculateBoundariesSurface();

    isMakeSuperStruct = true;   /* grid is fully built */
    return true;
}

