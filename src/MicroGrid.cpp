#include "MicroGrid.h"
#include "MathMatrix.h"
#include "MathGeomObjects.h"
#include <vector>
#include <array>
#include <map>
#include <algorithm>

ListMicroTetraedr::iterator JumpIterator(ListMicroTetraedr::iterator it, uint32_t count)
{
    uint32_t i;
    if(it != EmptyIterator())
        for(i=0; i<count; i++) it++;
    return it;
}

ListMicroTetraedr::iterator JumpIteratorBack(ListMicroTetraedr::iterator it, uint32_t count)
{
    uint32_t i;
    if(it != EmptyIterator())
        for(i=0; i<count; i++) it--;
    return it;
}

std::list<MicroNode>::iterator JumpIterator(std::list<MicroNode>::iterator it, uint32_t count)
{
    uint32_t i;
    for(i=0; i<count; i++) it++;
    return it;
}

std::list<MicroNode>::iterator JumpIteratorBack(std::list<MicroNode>::iterator it, uint32_t count)
{
    uint32_t i;
    for(i=0; i<count; i++) it--;
    return it;
}


MicroGrid::MicroGrid()
{
    setBeginPoint(MathPoint3D(0.0, 0.0, 0.0));
    setSizeGrid(0.0, 0.0, 0.0);
    isClearTetraedrs = false;
    isMakeSuperStruct = false;
    isLinkNodes = true;
    isNumNodesTetraedrs = false;
}


MicroGrid::MicroGrid(const MathPoint3D &begin, const double &dx, const double &dy, const double &dz)
{
    setBeginPoint(begin);
    setSizeGrid(dx, dy, dz);
    isClearTetraedrs = false;
    isMakeSuperStruct = false;
    isLinkNodes = false;
    isNumNodesTetraedrs = false;
}

void MicroGrid::setSizeGrid(const double &dx, const double &dy, const double &dz)
{
    superStructDx = dx;
    superStructDy = dy;
    superStructDz = dz;
    clear();
}

double MicroGrid::getSizeX()
{
    return superStructDx;
}

double MicroGrid::getSizeY()
{
    return superStructDy;
}

double MicroGrid::getSizeZ()
{
    return superStructDz;
}

void MicroGrid::setBeginPoint(const MathPoint3D &point)
{
    beginPoint = point;
    clear();
}

MathPoint3D MicroGrid::getBeginPoint()
{
    return beginPoint;
}

MathPoint3D MicroGrid::getBeginModelPoint()
{
    return beginModelPoint;
}

MathPoint3D MicroGrid::getEndModelPoint()
{
    return endModelPoint;
}

void MicroGrid::makeSuperStruct()
{
    std::list<MicroNode>::iterator itNode;
    uint32_t i;
    double constDx = beginPoint.getX() + superStructDx/2.0;
    double constDy = beginPoint.getY() + superStructDy/2.0;
    double constDz = beginPoint.getZ() + superStructDz/2.0;

    clear();
    for(i=0; i<4; i++) listNodes.push_front(MicroNode());
    itNode = listNodes.begin();
    itNode->point() = MathPoint3D(constDx, constDy - KOEF_SCALE_TETRAEDR*superStructDy/2.0, constDz + 2*KOEF_SCALE_TETRAEDR*superStructDz); itNode++;
    itNode->point() = MathPoint3D(constDx - KOEF_SCALE_TETRAEDR*superStructDx - 2.0*KOEF_SCALE_TETRAEDR*superStructDy, constDy - KOEF_SCALE_TETRAEDR*superStructDy/2.0, constDz - KOEF_SCALE_TETRAEDR*superStructDz/2.0); itNode++;
    itNode->point() = MathPoint3D(constDx + KOEF_SCALE_TETRAEDR*superStructDx + 2.0*KOEF_SCALE_TETRAEDR*superStructDy, constDy - KOEF_SCALE_TETRAEDR*superStructDy/2.0, constDz - KOEF_SCALE_TETRAEDR*superStructDz/2.0); itNode++;
    itNode->point() = MathPoint3D(constDx, constDy + KOEF_SCALE_TETRAEDR*superStructDy + KOEF_SCALE_TETRAEDR*superStructDx, constDz - KOEF_SCALE_TETRAEDR*superStructDz/2.0); itNode++;

    listTetraedr.push_front(MicroTetraedr());
    itNode = listNodes.begin();
    for(i=0; i<4; i++)
    {
        listTetraedr.begin()->nodes(i) = itNode;
        itNode->addFlags(NODE_IS_BOUNDARY | NODE_SUPERSTRUCT);
        itNode++;
    }
    listTetraedr.begin()->makeCircumSphere();
    listTetraedr.begin()->calculateBoundariesSurface();

    isMakeSuperStruct = true;
    isClearTetraedrs = false;
    isNumNodesTetraedrs = false;
}

bool MicroGrid::deleteSuperStruct()
{
    uint32_t i;
    ListMicroTetraedr::iterator itTetraedr;
    bool isSuperStructTetraedr;

    if(isClearTetraedrs) return false;
    if(isLinkNodes) return false;
    if(!isMakeSuperStruct) return false;
    isNumNodesTetraedrs = false;
    itTetraedr=listTetraedr.begin();
    while(itTetraedr!=listTetraedr.end())
    {
        isSuperStructTetraedr = false;
        for(i=0; i<4; i++)
            if(itTetraedr->nodes(i)->isFlags(NODE_SUPERSTRUCT))
            {
                isSuperStructTetraedr = true;
                break;
            }
        if(isSuperStructTetraedr)
        {
            for(i=0; i<4; i++)
                if(!itTetraedr->nodes(i)->isFlags(NODE_SUPERSTRUCT))
                    itTetraedr->nodes(i)->addFlags(NODE_IS_BOUNDARY);
            deleteTetraedr(itTetraedr);
        }
        else itTetraedr++;
    }
    isMakeSuperStruct = false;
    return true;
}

bool MicroGrid::addNode(const MicroNode &node)
{
    ListMicroTetraedr::iterator itTetraedr;
    std::list<ListMicroTetraedr::iterator> listDeleteTetraedrs;
    std::list<ListMicroTetraedr::iterator>::iterator listDeleteIterator, listDeleteIterator2;
    std::list<MicroSurface> listCreateTetraedrs;
    std::list<MicroSurface>::iterator listCreateIterator;
    std::vector<ListMicroTetraedr::iterator> vectorDoneTetraedrs;
    std::list<MicroNode>::iterator pushNode;
    std::list<MicroTriangle>::iterator itTriangle;
    uint32_t i, j, k;
    uint32_t numSharedSurface;
    bool deleteOk;
    bool createOk;

    /* Если точка не попадает в имеющуюся область */
    if(((MicroNode)node).point().getX() > beginPoint.getX() + superStructDx ||
            ((MicroNode)node).point().getY() > beginPoint.getY() + superStructDy ||
            ((MicroNode)node).point().getZ() > beginPoint.getZ() + superStructDz) return false;
    if(((MicroNode)node).point().getX() < beginPoint.getX() ||
            ((MicroNode)node).point().getY() < beginPoint.getY() ||
            ((MicroNode)node).point().getZ() < beginPoint.getZ()) return false;
    if(!listTetraedr.findTetraedr(((MicroNode)node).point(), itTetraedr)) return false;       /* невозможно добавить точку */

    for(i=0; i<4; i++)                                          /* если новый узел находится слишком близко к другим его не вставлять в структуру */
        if(((MicroNode)node).point() == itTetraedr->nodes(i)->point())
            return false;
    isNumNodesTetraedrs = false;
    listNodes.push_back(node);
    pushNode = listNodes.end();
    pushNode--;
    /* Обновление данных о размере модели */
    if(listNodes.size() > 5)
    {
        if(pushNode->point().getX() < beginModelPoint.getX()) beginModelPoint.setX(pushNode->point().getX());
        if(pushNode->point().getY() < beginModelPoint.getY()) beginModelPoint.setY(pushNode->point().getY());
        if(pushNode->point().getZ() < beginModelPoint.getZ()) beginModelPoint.setZ(pushNode->point().getZ());
        if(pushNode->point().getX() > endModelPoint.getX()) endModelPoint.setX(pushNode->point().getX());
        if(pushNode->point().getY() > endModelPoint.getY()) endModelPoint.setY(pushNode->point().getY());
        if(pushNode->point().getZ() > endModelPoint.getZ()) endModelPoint.setZ(pushNode->point().getZ());
    }
    else
    {
        beginModelPoint = pushNode->point();
        endModelPoint = pushNode->point();
    }
    /* этап удаления тетраэдров и образования полости в структуре */
    listDeleteTetraedrs.push_back(itTetraedr);
    listDeleteIterator = listDeleteTetraedrs.begin();
    while(listDeleteIterator != listDeleteTetraedrs.end())      /* пока список тетраэдров на удаление не станет пустым */
    {
        itTetraedr = *listDeleteIterator;
        for(i=0; i<4; i++)                                      /* проход по 4-м соседним тетраэдрам */
        {
            if(itTetraedr->tetraedrs(i) == EmptyIterator())     /* если нет одного из соседних тетраэдров */
            {
                if(itTetraedr->isBoundarySurface(i))           /* если itTetraedr имеет грань на границе структуры */
                    listTriangle.push_back(MicroTriangle(itTetraedr->nodes(i), itTetraedr->nodes(i + 1), itTetraedr->nodes(i + 2)));
                continue;
            }
            /* поиск номера плоскости у соседнего тетраэдра, граничащего с itTetraedr */
            if(!itTetraedr->tetraedrs(i)->getNumSurface(itTetraedr, numSharedSurface))
            {
                printf("Error in addNode: not found iterator (delete)\n");
                exit(-1);
            }
            if(itTetraedr->tetraedrs(i)->isInSphere(pushNode->point()) ||                               /* если соседний тетраэдр нарушает условие Делоне */
                    IsSingularTetraedr(pushNode->point(),                                               /* если тетраэдр, смежный с соседним будет вырожденным */
                                       itTetraedr->tetraedrs(i)->nodes(numSharedSurface)->point(),
                                       itTetraedr->tetraedrs(i)->nodes(numSharedSurface + 1)->point(),
                                       itTetraedr->tetraedrs(i)->nodes(numSharedSurface + 2)->point(),
                                       MIN_VOLUME_TETRAEDR) ||
                    itTetraedr->tetraedrs(i)->isInTetraedr(pushNode->point(), numSharedSurface) )
            {
                /* тетраэдр подлежит удалению */
                deleteOk = true;
                for(listDeleteIterator2=listDeleteTetraedrs.begin(); listDeleteIterator2!=listDeleteTetraedrs.end(); listDeleteIterator2++)
                    if(*listDeleteIterator2 == itTetraedr->tetraedrs(i))
                    {
                        /* если тетраэдр уже добавлен для удаления */
                        deleteOk = false;
                        break;
                    }
                if(deleteOk)
                {
                    /* удаление тетраэдра из списка для построения новых тетраэдров */
                    listCreateIterator = listCreateTetraedrs.begin();
                    while(listCreateIterator != listCreateTetraedrs.end())
                    {
                        if(listCreateIterator->tetraedr() == itTetraedr->tetraedrs(i))
                            listCreateIterator = listCreateTetraedrs.erase(listCreateIterator);
                        else listCreateIterator++;
                    }
                    listDeleteTetraedrs.push_back(itTetraedr->tetraedrs(i));
                }
            }
            else
            {
                /* грань тетраэдра используется для построения нового тетраэдра */
                createOk = true;
#ifdef DEBUG_MODE
                for(listCreateIterator=listCreateTetraedrs.begin(); listCreateIterator!=listCreateTetraedrs.end(); listCreateIterator++)
                    if(listCreateIterator->tetraedr() == itTetraedr->tetraedrs(i) && listCreateIterator->numSurface() == numSharedSurface)
                    {
                        printf("Error in addNode: 2 equal elements in listCreateTetraedrs\n");
                        exit(-1);
                    }
#endif
                for(listDeleteIterator2=listDeleteTetraedrs.begin(); listDeleteIterator2!=listDeleteTetraedrs.end(); listDeleteIterator2++)
                    if(*listDeleteIterator2 == itTetraedr->tetraedrs(i))
                    {
                        /* если тетраэдр уже помечен для удаления */
                        createOk = false;
                        break;
                    }
                if(createOk) listCreateTetraedrs.push_back(MicroSurface(itTetraedr->tetraedrs(i), numSharedSurface));
            }
        }
        listDeleteIterator = listDeleteTetraedrs.erase(listDeleteIterator);
        deleteTetraedr(itTetraedr);
    }
    /* Создание новых тетраэдров в полости, используя как грани другие тетраэдры */
    listCreateIterator = listCreateTetraedrs.begin();
    while(listCreateIterator != listCreateTetraedrs.end())
    {
        listTetraedr.push_front(MicroTetraedr());
        itTetraedr = listTetraedr.begin();

        itTetraedr->nodes(0) = pushNode;
        for(j=0; j<3; j++) itTetraedr->nodes(j + 1) = listCreateIterator->tetraedr()->nodes(listCreateIterator->numSurface() + j);
        for(j=0; j<4; j++)
        {
            if(j == 1)
            {
                itTetraedr->tetraedrs(j) = listCreateIterator->tetraedr();
                if(!listCreateIterator->tetraedr()->getNumSurface(itTetraedr->nodes(j)->point(),
                                                    itTetraedr->nodes(j + 1)->point(),
                                                    itTetraedr->nodes(j + 2)->point(), numSharedSurface))
                {
                    printf("Error in addNode: not found iterator (add)\n");
                    exit(-1);
                }
                listCreateIterator->tetraedr()->tetraedrs(numSharedSurface) = itTetraedr;
            }
            else
            {
                for(k=0; k<vectorDoneTetraedrs.size(); k++)
                    if(vectorDoneTetraedrs[k]->isNode3Tetraedr(itTetraedr->nodes(j)->point(),
                                                               itTetraedr->nodes(j + 1)->point(),
                                                               itTetraedr->nodes(j + 2)->point()))
                    {
                        itTetraedr->tetraedrs(j) = vectorDoneTetraedrs[k];
                        if(!vectorDoneTetraedrs[k]->getNumSurface(itTetraedr->nodes(j)->point(),
                                                              itTetraedr->nodes(j + 1)->point(),
                                                              itTetraedr->nodes(j + 2)->point(), numSharedSurface))
                        {
                            printf("Error in addNode: not found iterator (add 2)\n");
                            exit(-1);
                        }
                        if(vectorDoneTetraedrs[k]->tetraedrs(numSharedSurface) != EmptyIterator())
                        {
                            printf("Error in addNode: not add new tetraedr\n");
                            exit(-1);
                        }
                        vectorDoneTetraedrs[k]->tetraedrs(numSharedSurface) = itTetraedr;
                        break;
                    }
            }
        }
        vectorDoneTetraedrs.push_back(itTetraedr);
        listCreateIterator = listCreateTetraedrs.erase(listCreateIterator);
    }
    /* Создание новых тетраэдров в полости, используя как грани треугольные плоскости на границе */
    itTriangle = listTriangle.begin();
    while(itTriangle != listTriangle.end())
    {
        listTetraedr.push_front(MicroTetraedr());
        itTetraedr = listTetraedr.begin();

        itTetraedr->nodes(0) = pushNode;
        for(j=0; j<3; j++) itTetraedr->nodes(j + 1) = itTriangle->node(j);
        if(!IsSingularTetraedr(itTetraedr->nodes(0)->point(),
                               itTetraedr->nodes(1)->point(),
                               itTetraedr->nodes(2)->point(),
                               itTetraedr->nodes(3)->point(),
                               MIN_VOLUME_TETRAEDR))
        {
            for(j=0; j<4; j++)
            {
                if(j == 1)
                {
                    itTetraedr->tetraedrs(j) = EmptyIterator();
                }
                else
                {
                    for(k=0; k<vectorDoneTetraedrs.size(); k++)
                        if(vectorDoneTetraedrs[k]->isNode3Tetraedr(itTetraedr->nodes(j)->point(),
                                                                   itTetraedr->nodes(j + 1)->point(),
                                                                   itTetraedr->nodes(j + 2)->point()))
                        {
                            itTetraedr->tetraedrs(j) = vectorDoneTetraedrs[k];
                            if(!vectorDoneTetraedrs[k]->getNumSurface(itTetraedr->nodes(j)->point(),
                                                                  itTetraedr->nodes(j + 1)->point(),
                                                                  itTetraedr->nodes(j + 2)->point(), numSharedSurface))
                            {
                                printf("Error in addNode: not found iterator (add 3)\n");
                                exit(-1);
                            }
                            vectorDoneTetraedrs[k]->tetraedrs(numSharedSurface) = itTetraedr;
                            break;
                        }
                }
            }
            vectorDoneTetraedrs.push_back(itTetraedr);
        }
        else
        {
            listTetraedr.erase(itTetraedr);
            pushNode->addFlags(NODE_IS_BOUNDARY);
        }
        itTriangle = listTriangle.erase(itTriangle);
    }
    for(i=0; i<vectorDoneTetraedrs.size(); i++)
    {
        vectorDoneTetraedrs[i]->makeCircumSphere();
        vectorDoneTetraedrs[i]->calculateBoundariesSurface();
    }

    return true;
}


bool MicroGrid::deleteOneTetraedr(ListMicroTetraedr::iterator &it)
{
    uint32_t i;

    if(isClearTetraedrs) return false;
    if(isLinkNodes) clearLinkNodes();
    isNumNodesTetraedrs = false;
    for(i=0; i<4; i++) it->nodes(i)->addFlags(NODE_IS_BOUNDARY);
    deleteTetraedr(it);
    return true;
}

bool MicroGrid::linkNodes()
{
    ListMicroTetraedr::iterator itTetraedr;
    std::list<MicroNode>::iterator itNode;
    uint32_t i, j;

    if(isLinkNodes) return true;
    if(isClearTetraedrs) return false;
    isLinkNodes = true;
    isNumNodesTetraedrs = false;
    /* связывание узлов */
    for(itTetraedr=listTetraedr.begin(); itTetraedr!=listTetraedr.end(); itTetraedr++)
        for(i=0; i<4; i++)
            for(j=0; j<4; j++)
            {
                if(i == j) continue;
                itTetraedr->nodes(i)->addNeighbourNode(itTetraedr->nodes(j));
            }
    /* удаление свободных узлов */
    itNode = listNodes.begin();
    while(itNode != listNodes.end())
    {
        if(!itNode->numNeighbourNodes())
        {
            itNode = listNodes.erase(itNode);
        }
        else itNode++;
    }
    return true;
}

void MicroGrid::clearLinkNodes()
{
    std::list<MicroNode>::iterator itNode;

    if(isLinkNodes) return;
    if(isClearTetraedrs) return;

    for(itNode = listNodes.begin(); itNode!=listNodes.end(); itNode++)
        itNode->clearNeighbourNodes();
    isLinkNodes = false;
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

bool MicroGrid::isAllDelone()
{
    uint32_t i;
    ListMicroTetraedr::iterator itTetraedr;
    bool isDelone;

    if(isClearTetraedrs) return false;
    if(!isMakeSuperStruct) return false;

    /* Проверка на условие Делоне */
    isDelone = true;
    for(itTetraedr = listTetraedr.begin(); itTetraedr != listTetraedr.end(); itTetraedr++)
    {
        for(i=0; i<4; i++)
        {
            if(itTetraedr->tetraedrs(i) == EmptyIterator()) continue;
            if(!IsDeloneTetraedrs(itTetraedr, itTetraedr->tetraedrs(i)))
            {
                isDelone = false;
                break;
            }
        }
    }
    return isDelone;
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

double MicroGrid::volumeSuperStruct()
{
    ListMicroTetraedr::iterator itTetraedr;
    double volume = 0.0;

    if(isClearTetraedrs) return 0.0;

    for(itTetraedr = listTetraedr.begin(); itTetraedr != listTetraedr.end(); itTetraedr++)
        volume += itTetraedr->volume();
    return volume;
}

void MicroGrid::clearTetraedrs()
{
    isClearTetraedrs = true;
    listTetraedr.clear();
}

void MicroGrid::clear()
{
    isClearTetraedrs = false;
    isMakeSuperStruct = false;
    isLinkNodes = false;
    isNumNodesTetraedrs = false;
    beginModelPoint = MathPoint3D(0.0, 0.0, 0.0);
    endModelPoint = MathPoint3D(0.0, 0.0, 0.0);
    listNodes.clear();
    listTetraedr.clear();
    listTriangle.clear();
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
        node.setPoint(MathPoint3D(mesh.nodes[i][0], mesh.nodes[i][1], mesh.nodes[i][2]));
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
        MathPoint3D p0 = nodeIt[t[0]]->point(), p1 = nodeIt[t[1]]->point();
        MathPoint3D p2 = nodeIt[t[2]]->point(), p3 = nodeIt[t[3]]->point();
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
        MathPoint3D p0 = nodeIt[t[0]]->point(), p1 = nodeIt[t[1]]->point();
        MathPoint3D p2 = nodeIt[t[2]]->point(), p3 = nodeIt[t[3]]->point();
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

    /* 5. Neighbour links (for GUI drawing) and sequential numbering. */
    linkNodes();
    setNumberNodesTetraedrs();

    /* 6. Mark boundary surfaces (needs NODE_IS_BOUNDARY set above). */
    for(ListMicroTetraedr::iterator it = listTetraedr.begin(); it != listTetraedr.end(); ++it)
        it->calculateBoundariesSurface();

    isMakeSuperStruct = true;   /* grid is fully built */
    return true;
}

bool MicroGrid::tradeTetraedrs(ListMicroTetraedr::iterator tetraedr1, ListMicroTetraedr::iterator tetraedr2, ListMicroTetraedr::iterator *newIterator)
{
    ListMicroTetraedr::iterator tableNeighbourTetraedrs[3];
    ListMicroTetraedr::iterator tetraedrs[3];
    MicroTetraedr oldTetraedrs[2];
    uint32_t i, j, k;
    uint32_t numFreeVertex1, numFreeVertex2;

    for(i=0; i<4; i++)
        if(!tetraedr2->isNodeTetraedr(tetraedr1->nodes(i)->point()))
        {
            numFreeVertex1 = i;
            break;
        }
#ifdef DEBUG_MODE
    if(i == 4)
    {
        printf("Error in trade tetraedrs \n");
        exit(-1);
    }
#endif

    for(i=0; i<4; i++)
        if(!tetraedr1->isNodeTetraedr(tetraedr2->nodes(i)->point()))
        {
            numFreeVertex2 = i;
            break;
        }
#ifdef DEBUG_MODE
    if(i == 4)
    {
        printf("Error in trade tetraedrs \n");
        exit(-1);
    }
#endif

    for(i=0; i<3; i++)
    {
#ifdef DEBUG_MODE
        /* проверка двух тетраэдров, что они смежные */
        for(j=0; j<3; j++) if(tetraedr1->nodes(numFreeVertex1 + ((i + 1) % 3) + 1)->point() == tetraedr2->nodes(numFreeVertex2 + j + 1)->point()) break;
        if(j == 3)
        {
            printf("Error in trade tetraedrs \n");
            exit(-1);
        }
        for(j=0; j<3; j++) if(tetraedr1->nodes(numFreeVertex1 + ((i + 2) % 3) + 1)->point() == tetraedr2->nodes(numFreeVertex2 + j + 1)->point()) break;
        if(j == 3)
        {
            printf("Error in trade tetraedrs \n");
            exit(-1);
        }
#endif
        /* создание таблицы соседних тетраэдров (внешних к создаваемым) из tetraedr2 */
        for(j=0; j<4; j++)
            if(tetraedr2->tetraedrs(j) != EmptyIterator())
                if(tetraedr2->tetraedrs(j)->isNode3Tetraedr(tetraedr1->nodes(numFreeVertex1 + ((i + 1) % 3) + 1)->point(),
                                                            tetraedr1->nodes(numFreeVertex1 + ((i + 2) % 3) + 1)->point(),
                                                            tetraedr2->nodes(numFreeVertex2)->point()))
                {
                    tableNeighbourTetraedrs[i] = tetraedr2->tetraedrs(j);
                    break;
                }

        if(j == 4) tableNeighbourTetraedrs[i] = EmptyIterator();
    }

    /* если трейд невозможно выполнить в принципе - то возврат false */
    for(i=0; i<4; i++)
        if(tetraedr1->tetraedrs(i) != EmptyIterator())
            if(tetraedr1->tetraedrs(i)->isNodeTetraedr(tetraedr1->nodes(numFreeVertex1)->point()) &&
                    tetraedr1->tetraedrs(i)->isNodeTetraedr(tetraedr2->nodes(numFreeVertex2)->point()))
                return false;
    for(i=0; i<4; i++)
        if(tetraedr2->tetraedrs(i) != EmptyIterator())
            if(tetraedr2->tetraedrs(i)->isNodeTetraedr(tetraedr2->nodes(numFreeVertex2)->point()) &&
                    tetraedr2->tetraedrs(i)->isNodeTetraedr(tetraedr1->nodes(numFreeVertex1)->point()))
                return false;

    /* добавление третьего тетраэдра в список */
    listTetraedr.push_back(MicroTetraedr());
    tetraedrs[0] = tetraedr1;
    tetraedrs[1] = tetraedr2;
    tetraedrs[2] = listTetraedr.end();
    tetraedrs[2]--;

    oldTetraedrs[0] = *tetraedr1;
    oldTetraedrs[1] = *tetraedr2;

    for(i=0; i<3; i++)
    {
        tetraedrs[i]->nodes(0) = oldTetraedrs[0].nodes(numFreeVertex1);
        tetraedrs[i]->nodes(1) = oldTetraedrs[0].nodes(numFreeVertex1 + ((i + 1) % 3) + 1);
        tetraedrs[i]->nodes(2) = oldTetraedrs[0].nodes(numFreeVertex1 + ((i + 2) % 3) + 1);
        tetraedrs[i]->nodes(3) = oldTetraedrs[1].nodes(numFreeVertex2);
    }

    /* обновление информации о соседних тетраэдрах */
    for(i=0; i<3; i++)
    {
        tetraedrs[i]->tetraedrs(0) = oldTetraedrs[0].tetraedrs(numFreeVertex1 + ((i + 2) % 4));
        tetraedrs[i]->tetraedrs(1) = tableNeighbourTetraedrs[i];
        tetraedrs[i]->tetraedrs(2) = tetraedrs[(i + 1) % 3];
        tetraedrs[i]->tetraedrs(3) = tetraedrs[(i + 2) % 3];
        for(j=0; j<2; j++)
            if(tetraedrs[i]->tetraedrs(j) != EmptyIterator())
                for(k=0; k<4; k++)    /* перебор по итераторам на соседние тетраэдры у соседа(не подвергавшегося трейду) */
                    if(tetraedrs[i]->tetraedrs(j)->tetraedrs(k) == tetraedr1 || tetraedrs[i]->tetraedrs(j)->tetraedrs(k) == tetraedr2)
                        tetraedrs[i]->tetraedrs(j)->tetraedrs(k) = tetraedrs[i];
    }

    /* вычисление описанной окружности */
    for(i=0; i<3; i++)
    {
        tetraedrs[i]->makeCircumSphere();
        tetraedrs[i]->calculateBoundariesSurface();
    }

    if(newIterator) *newIterator = tetraedrs[2];
    /* трейд выполнен успешно */
    return true;
}


bool IsObjectVector(std::vector<std::list<MicroNode>::iterator> &list, std::list<MicroNode>::iterator &it)
{
    uint32_t i;
    for(i=0; i<list.size(); i++)
        if(list[i] == it) return true;
    return false;
}


/* Функция удаления тетраэдра из сетки без остатков */
void MicroGrid::deleteTetraedr(ListMicroTetraedr::iterator &it)
{
    uint32_t i, numSurface;

    /* Обнуление итераторов на данный тетраэдр у соседей */
    for(i=0; i<4; i++)
    {
        if(it->tetraedrs(i) != EmptyIterator())
        {
            if(!it->tetraedrs(i)->getNumSurface(it, numSurface))
            {
                printf("Error in deleteTetraedr: not found iterator\n");
                exit(-1);
            }
            it->tetraedrs(i)->tetraedrs(numSurface) = EmptyIterator();
        }
    }
    it = listTetraedr.erase(it);
}

//EOF
