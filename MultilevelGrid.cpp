/***************************************************************************
    File        : MultilevelGrid.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Jan 18 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "MultilevelGrid.h"
#include "Hexahedron.h"
using namespace std;
namespace grd {


// Description
//
MultilevelGrid::MultilevelGrid()
{
  topLevel     = 0;
  bndTopLevel  = 0;
  noOfFaceSets = 0;

  // create grid level 0
  GridLevel* grid0 = new GridLevel;

  // put grid level 0 into grid stack
  grid.push(grid0);
  noOfVertices.push(grid0->noOfVertices());
  noOfEdges.push(grid0->noOfEdges());
  noOfTriangles.push(grid0->noOfTriangles());
  noOfQuadrangles.push(grid0->noOfQuadrangles());
  noOfTetrahedra.push(grid0->noOfTetrahedra());
  noOfOctahedra.push(grid0->noOfOctahedra());
  noOfHexahedra.push(grid0->noOfHexahedra());
  noOfPrisms.push(grid0->noOfPrisms());
  noOfPyramids.push(grid0->noOfPyramids());

}


// Description
//
MultilevelGrid::~MultilevelGrid()
{
  int noOfGrids = grid.size();
  for (int i = 0; i < noOfGrids; i++)
    delete grid[i];
}



// Description
// same algorithm as refine volume
// process only boundary elements
void
MultilevelGrid::refineBoundary()
{
  int i;

  // check if next grid level already exits
  if (topLevel <= bndTopLevel) {
    GridLevel* newGridLevel = new GridLevel;
    grid.push(newGridLevel);
  }

  // top-down
  // refine marked elements
  // irregular elements are refined regular
  for (i = bndTopLevel; i >= 0; i--) {
    grid[i]->refineBoundary(*grid[i+1]);
  }

  // remove non used vertices and edges
  for (i = 0; i <= (bndTopLevel+1); i++)
    grid[i]->refVerticesAndEdges();

  for (i = 0; i <= (bndTopLevel+1); i++)
    grid[i]->removeUnrefVerticesAndEdges();

  for (i = 0; i <= (bndTopLevel+1); i++)
    grid[i]->unrefVerticesAndEdges();

  //bottom-up
  // coarsen grid level
  for(i = 1; i <= bndTopLevel; i++) {
    grid[i]->coarsenBoundary();
  }
  //close grid level
  for(i = 0; i < bndTopLevel; i++) {
    grid[i]->closeBoundary();
  }

  // compute no. of levels
  if (topLevel <= bndTopLevel) {
    if (!grid[bndTopLevel+1]->empty()) {
      bndTopLevel++;
      noOfVertices.push(0);
      noOfEdges.push(0);
      noOfTriangles.push(0);
      noOfQuadrangles.push(0);
    }
    else if (bndTopLevel > 0) {
      delete grid[bndTopLevel+1]; // remove bndTopLevle + 1
      grid.pop();
      noOfVertices.pop();
      noOfEdges.pop();
      noOfTriangles.pop();
      noOfQuadrangles.pop();
      if (grid[bndTopLevel]->empty()) {
        delete (grid[bndTopLevel]); // remove bndTopLevel
        grid.pop();
        noOfVertices.pop();
        noOfEdges.pop();
        noOfTriangles.pop();
        noOfQuadrangles.pop();
        bndTopLevel--;
      }
    }
  }
  else
    bndTopLevel++;

  // recompute number of elements
  updateNoOfElements();
}



// Description
// top-down / bottom-up algorithm
void
MultilevelGrid::refine()
{
  int i;

  // Create new grid level and push into stack
  if (bndTopLevel <= topLevel) {
    GridLevel* newGridLevel = new GridLevel;
    grid.push(newGridLevel);
  }

  // top-down
  // refine marked elements
  // irregular elements are refined regular
  for (i = topLevel; i >= 0; i--) {
    grid[i]->refine(*grid[i+1]);
  }

  // remove unused vertices and edges from lists
  for (i = 0; i <= (topLevel+1); i++)
    grid[i]->refVerticesAndEdges();

  for (i = 0; i <= (topLevel+1); i++)
    grid[i]->removeUnrefVerticesAndEdges();

  for (i = 0; i <= (topLevel+1); i++)
    grid[i]->unrefVerticesAndEdges();

  //bottom-up
  // coarsen grid level
  for(i = 1; i <= topLevel; i++) {
    grid[i]->coarsen();
  }
  //close grid level
  for(i = 0; i <= topLevel; i++) {
    grid[i]->close(*grid[i+1]);
  }

  // compute new no. of levels
  if (bndTopLevel <= topLevel)
  {
    if (!grid[topLevel+1]->empty())
    {
      noOfVertices.push(0);
      noOfEdges.push(0);
      noOfTriangles.push(0);
      noOfQuadrangles.push(0);
      noOfTetrahedra.push(0);
      noOfOctahedra.push(0);
      noOfHexahedra.push(0);
      noOfPrisms.push(0);
      noOfPyramids.push(0);
      topLevel++;
    }
    else if (topLevel > 0) {
      delete grid[topLevel+1];
      grid.pop();
      noOfVertices.pop();
      noOfEdges.pop();
      noOfTriangles.pop();
      noOfQuadrangles.pop();
      noOfTetrahedra.pop();
      noOfOctahedra.pop();
      noOfHexahedra.pop();
      noOfPrisms.pop();
      noOfPyramids.pop();
      if (grid[topLevel]->empty()) {
        delete grid[topLevel];
        grid.pop();
        noOfVertices.pop();
        noOfEdges.pop();
        noOfTriangles.pop();
        noOfQuadrangles.pop();
        noOfTetrahedra.pop();
        noOfOctahedra.pop();
        noOfHexahedra.pop();
        noOfPrisms.pop();
        noOfPyramids.pop();
        topLevel--;
      }
    }
  }
  else
    topLevel++;
  // Recompute number of elements
  updateNoOfElements();


  for (i = 0; i <= topLevel; i++)
  {
    list<Vertex*>::iterator p;
    list<Vertex*>* lt = grid[i]->getVertexList();
    for (p = lt->begin(); p != lt->end(); ++p)
      (*p)->setId(-9);
  }

  for (i = 0; i <= topLevel; i++)
  {
    list<Edge*>::iterator p;
    list<Edge*>* lt = grid[i]->getEdgeList();
    for (p = lt->begin(); p != lt->end(); ++p)
    {
      for (int vv=0; vv < 2; vv++)
      {
        Vertex* vertex = (*p)->getVertex(vv);
        int vid = vertex->getId();
        if (vid != -9)
          cerr << "ERROR: wrong id: " << vid << "  for level: " << i << '\n';
      }
    }
  }
  for (i = 0; i <= topLevel; i++)
  {
    list<Element*>::iterator p;
    list<Element*>* lt = grid[i]->getHexahedronList();
    for (p = lt->begin(); p != lt->end(); ++p)
    {
      if ((*p)->isIrregular())
      {
        Vertex* c = ((Hexahedron*)(*p))->getBarycenter();
        if (c)
        {
          int vid = c->getId();
          if (vid != -9)
            cerr << "ERROR: barycenter wrong id: " << vid << "  for level: " << i << '\n';
        }
        else
          cerr << "ERROR: c == 0\n";
      }
    }
  }

  // set Grid level
  bndTopLevel = topLevel;
}





//
// Description
// compute mesh topology from a list of vertices and a list of elements
void
MultilevelGrid::buildCoarseMesh(vector<Vertex*>& vertex,vector<Element*>& element)
{
  int i,j,k,l,m;
  int face;
  int edge;
  int setCounter;
  int noOfVertices = vertex.size();
  int noOfElements = element.size();
  Vertex* vt;
  Vertex* vt1;
  Vertex* vt2;
  Vertex* nvt1;
  Vertex* nvt2;
  Edge*   theEdge;
  Edge*   nedge;
  Element* theElement;
  Element* neigh1;
  Element* neigh2;

  list<Vertex*>*  vertexList = grid[0]->getVertexList();
  list<Edge*>*    edgeList   = grid[0]->getEdgeList();
  list<Element*>* triList    = grid[0]->getTriangleList();
  list<Element*>* quaList    = grid[0]->getQuadrangleList();
  list<Element*>* tetList    = grid[0]->getTetrahedronList();
  vector<Stack<Element*>*> vtList(noOfVertices);
  Stack<Element*>* stack1;
  Stack<Element*>* stack2;


  // set vertex and element list
  for (i = 0; i < noOfVertices; i++)
    vertexList->push_back(vertex[i]);
  for (i = 0; i < noOfElements; i++)
  {
    int type = element[i]->type();
    switch(type)
    {
      case GRD_TRIANGLE:
        triList->push_back(element[i]);
        break;
      case GRD_QUADRANGLE:
        quaList->push_back(element[i]);
        break;
      case GRD_TETRAHEDRON:
        tetList->push_back(element[i]);
        break;
      default:
        break;
    }
  }

  //===============================//
  // get neighborhood relations   //
  //==============================//

  // For each vertex search the elements
  // which share it, use an stack
  for (i = 0; i < noOfVertices; i++)
    vtList[i] = new Stack<Element*>;

  for (i= 0; i < noOfElements; i++) {
    for (j = 0; j < element[i]->getNoOfVertices(); j++) {
      vt = element[i]->getVertex(j);
      vtList[vt->getId() - 1]->push(element[i]);
    }
  }


  // Create element neighborhood relations
  for (l = 0; l < noOfElements; l++)
  {
    theElement = element[l];
    int type = theElement->type();
    if (type == GRD_TRIANGLE || type == GRD_QUADRANGLE)
      continue;

    for (face = 0; face < theElement->getNoOfFaces(); face++)
    {
      // for each vertex at face
      int noOfVertexAtFace = theElement->noOfVertexAtFace(face);
      Stack<Element*>** stck = new Stack<Element*>*[noOfVertexAtFace];
      for (m = 0; m < noOfVertexAtFace; m++)
      {
        vt = theElement->getVertexAtFace(face,m);
        stck[m] = vtList[vt->getId() - 1];
      }

      // Init:
      // intersection of first two lists
      list<Element*>* lst1 = new list<Element*>;
      list<Element*>* lst2 = new list<Element*>;
      for (i = 0; i < stck[0]->size(); i++)
      {
        neigh1 = (*stck[0])[i];
        if (neigh1 != theElement)
        {
          for (j = 0; j < stck[1]->size(); j++)
          {
            neigh2 = (*stck[1])[j];
            if (neigh2 == neigh1)
            {
              lst1->push_back(neigh1);
            }
          }
        }
      }
      // Intersect with restant lists
      for (m = 2; m < noOfVertexAtFace; m++)
      {
        // intersect allways two lists
        for (i = 0; i < stck[m]->size(); i++)
        {
          neigh1 = (*stck[m])[i];
          typedef list<Element*>::iterator ElmI;
          for (ElmI p = lst1->begin(); p != lst1->end(); p++)
          {
            neigh2 = *p;
            if (neigh1 == neigh2)
            {
              lst2->push_back(neigh1);
            }
          }

        }
        lst1->clear();
        list<Element*>* tmpLst = lst1;
        lst1 = lst2;
        lst2 = tmpLst;
      }

      // Set elements face
      if (lst1->size() == 1)
        theElement->setNeighbor(face,lst1->back());
      else if (lst1->size() > 1)
        cerr << "ERROR: more than one neighbor for a face\n";

      // remove auxiliary data
      delete lst1;
      delete lst2;
      delete [] stck;

    } // for each face
  } // for each element

  /* create edge list */
  setCounter = 0;
  int totNoOfEdges = 0;
  //  while((theElement = ++elIter1)) {
  for (l = 0; l < noOfElements; l++) {
    theElement = element[l];
    if (setCounter > 1) {
      cerr << "panic! set counter > 1\n";
      exit(1);
    }
    for (edge = 0; edge < theElement->getNoOfEdges(); edge++) {
      setCounter = 0;
      theEdge = theElement->getEdge(edge);
      if (!theEdge) {
        vt1 = theElement->getVertexAtEdge(edge,0);
        stack1 = vtList[vt1->getId() - 1];
        vt2 = theElement->getVertexAtEdge(edge,1);
        stack2 = vtList[vt2->getId() - 1];

        for (i = 0; i < stack1->size(); i++) {
          if (!(theElement->getEdge(edge))) {
            neigh1 = (*stack1)[i];
            if (neigh1 != theElement) {
              for (j = 0; j < stack2->size(); j++) {
                if (!(theElement->getEdge(edge))) {
                  neigh2 = (*stack2)[j];
                  if (neigh1 == neigh2) {
                    for (k = 0; k < neigh1->getNoOfEdges(); k++) {
                      nedge = neigh1->getEdge(k);
                      if (nedge) {
                        nvt1 = nedge->getVertex(0);
                        nvt2 = nedge->getVertex(1);
                        if (vt1 == nvt1 && vt2 == nvt2) {
                          theElement->setEdge(edge,nedge);
                          setCounter++;
                          break;
                        }
                        else if (vt1 == nvt2 && vt2 == nvt1) {
                          theElement->setEdge(edge,nedge);
                          setCounter++;
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (setCounter == 0) { //(!(elem->getEdge(edge)))
        theEdge = new Edge;
        edgeList->push_back(theEdge);
        vt1 = theElement->getVertexAtEdge(edge,0);
        vt2 = theElement->getVertexAtEdge(edge,1);
        theEdge->setVertices(vt1,vt2);
        theElement->setEdge(edge,theEdge);
        totNoOfEdges++;
      }
    }
  }

  // Extract the boundary
  grid[0]->extractBoundary();

  // Recompute the number of elements
  updateNoOfElements();
}



// Description
// Recompute the number of elements
void
MultilevelGrid::updateNoOfElements()
{
  int i;
  int counter;
  int upLevel = (topLevel < bndTopLevel) ? bndTopLevel : topLevel;

  // recompute number of vertices
  noOfVertices[0] = grid[0]->noOfVertices();
  for (i = 1; i <= upLevel; i++)
    noOfVertices[i]    = noOfVertices[i-1] + grid[i]->noOfVertices();

  // recompute no. of edges
  noOfEdges[0] = grid[0]->noOfEdges();
  typedef list<Edge*>::iterator EdgI;
  counter = 0;
  for (i = 0; i < upLevel; i++) {
    list<Edge*>* lt = grid[i]->getEdgeList();
    for (EdgI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!((*p)->isRefined()))
        counter++;
    }
    noOfEdges[i+1] = counter + grid[i+1]->noOfEdges();
  }

  // recompute no. elements
  noOfTriangles[0]   = grid[0]->noOfTriangles();
  noOfQuadrangles[0] = grid[0]->noOfQuadrangles();
  noOfTetrahedra[0]  = grid[0]->noOfTetrahedra();
  noOfOctahedra[0]   = grid[0]->noOfOctahedra();
  noOfHexahedra[0]   = grid[0]->noOfHexahedra();
  noOfPrisms[0]      = grid[0]->noOfPrisms();
  noOfPyramids[0]    = grid[0]->noOfPyramids();


  typedef list<Element*>::iterator EI;
  list<Element*>* lt;

  counter = 0;
  for (i = 0; i < upLevel; i++) {
    // Triangles
    lt = grid[i]->getTriangleList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfTriangles[i+1] = counter + grid[i+1]->noOfTriangles();
  }
  counter = 0;
  for (i = 0; i < upLevel; i++) {
    // Quadriangle
    lt = grid[i]->getQuadrangleList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfQuadrangles[i+1] = counter + grid[i+1]->noOfQuadrangles();
  }
  counter = 0;
  for (i = 0; i < topLevel; i++) {
    // Tetrahedron
    lt = grid[i]->getTetrahedronList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfTetrahedra[i+1] = counter + grid[i+1]->noOfTetrahedra();
  }
  counter = 0;
  for (i = 0; i < topLevel; i++) {
    // Octahedron
    lt = grid[i]->getOctahedronList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfOctahedra[i+1] = counter + grid[i+1]->noOfOctahedra();
  }
  counter = 0;
  for (i = 0; i < topLevel; i++) {
    // Hexahedron
    lt = grid[i]->getHexahedronList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfHexahedra[i+1] = counter + grid[i+1]->noOfHexahedra();
  }
  counter = 0;
  for (i = 0; i < topLevel; i++) {
    // Prism
    lt = grid[i]->getPrismList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfPrisms[i+1] = counter + grid[i+1]->noOfPrisms();
  }
  counter = 0;
  for (i = 0; i < topLevel; i++) {
    // Pyramid
    lt = grid[i]->getPyramidList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
    noOfPyramids[i+1] = counter + grid[i+1]->noOfPyramids();
  }
}


//
//
void
MultilevelGrid::buildMultilevelBoundary(Stack<double**>& mlVertex,Stack<int**>& mlElement)
{
  int counter = 0;
  list<Element*>* lt = grid[0]->getTriangleList();
  Element* tri;
  typedef list<Element*>::iterator TriI;
  double** pos = mlVertex[0];
  int**    elm = mlElement[0];

  int noOfLevel = mlVertex.size();
  list<TriSkeleton*> triSk1;
  list<TriSkeleton*> triSk2;

  for (TriI p = lt->begin(); p != lt->end(); ++p)
  {
    tri = *p;
    TriSkeleton* skl = new TriSkeleton;
    TriSkeleton* ch0 = new TriSkeleton;
    TriSkeleton* ch1 = new TriSkeleton;
    TriSkeleton* ch2 = new TriSkeleton;
    TriSkeleton* ch3 = new TriSkeleton;

    for (int n = 0; n < 3; n++)
    {
      ch0->pos[0][n] = pos[elm[counter + 0][0]][n];
      ch0->pos[1][n] = pos[elm[counter + 0][1]][n];
      ch0->pos[2][n] = pos[elm[counter + 0][2]][n];

      ch1->pos[0][n] = pos[elm[counter + 1][0]][n];
      ch1->pos[1][n] = pos[elm[counter + 1][1]][n];
      ch1->pos[2][n] = pos[elm[counter + 1][2]][n];

      ch2->pos[0][n] = pos[elm[counter + 2][0]][n];
      ch2->pos[1][n] = pos[elm[counter + 2][1]][n];
      ch2->pos[2][n] = pos[elm[counter + 2][2]][n];

      ch3->pos[0][n] = pos[elm[counter + 3][0]][n];
      ch3->pos[1][n] = pos[elm[counter + 3][1]][n];
      ch3->pos[2][n] = pos[elm[counter + 3][2]][n];
    }

    // set position of root TriSkeleton
    for (int i = 0; i < 3; i++)
    {
      Vertex* v = tri->getVertex(i);
      skl->pos[i][0] = (*v)[0];
      skl->pos[i][1] = (*v)[1];
      skl->pos[i][2] = (*v)[2];
    }
    // set children
    skl->child[0] = ch0;
    skl->child[1] = ch1;
    skl->child[2] = ch2;
    skl->child[3] = ch3;

    ((Triangle*)tri)->setSkeleton(skl);
    counter += 4;

    triSk1.push_back(ch0);
    triSk1.push_back(ch1);
    triSk1.push_back(ch2);
    triSk1.push_back(ch3);
  }

  list<TriSkeleton*>::iterator q;
  for (int level = 1; level < noOfLevel; level++)
  {
    cerr << "Processing level: " << level << '\n';
    cerr << "List size: " << (triSk1.size()) << '\n';
    pos = mlVertex[level];
    elm = mlElement[level];
    counter = 0;
    int chIdx = 0;
    for (q = triSk1.begin(); q != triSk1.end(); ++q)
    {
      TriSkeleton* skl = *q;

      TriSkeleton* ch0 = new TriSkeleton;
      TriSkeleton* ch1 = new TriSkeleton;
      TriSkeleton* ch2 = new TriSkeleton;
      TriSkeleton* ch3 = new TriSkeleton;

      int chLoc = (chIdx%4);
      if (chLoc == 3) chLoc = 0;
      int chCh = 0;
      for (int n = 0; n < 3; n++)
      {
        ch0->pos[0][n] = pos[elm[counter][0]][n];
        ch0->pos[1][n] = pos[elm[counter][1]][n];
        ch0->pos[2][n] = pos[elm[counter][2]][n];
        chCh++;
        ch1->pos[0][n] = pos[elm[counter + 1][0]][n];
        ch1->pos[1][n] = pos[elm[counter + 1][1]][n];
        ch1->pos[2][n] = pos[elm[counter + 1][2]][n];
        chCh++;
        ch2->pos[0][n] = pos[elm[counter + 2][0]][n];
        ch2->pos[1][n] = pos[elm[counter + 2][1]][n];
        ch2->pos[2][n] = pos[elm[counter + 2][2]][n];
        chCh++;
        ch3->pos[0][n] = pos[elm[counter + 3][0]][n];
        ch3->pos[1][n] = pos[elm[counter + 3][1]][n];
        ch3->pos[2][n] = pos[elm[counter + 3][2]][n];
      }

      // set children
      skl->child[0] = ch0;
      skl->child[1] = ch1;
      skl->child[2] = ch2;
      skl->child[3] = ch3;

      triSk2.push_back(ch0);
      triSk2.push_back(ch1);
      triSk2.push_back(ch2);
      triSk2.push_back(ch3);
      counter += 4;
      chIdx++;
    }

    // copy lists
    triSk1.clear();
    for (q = triSk2.begin(); q != triSk2.end(); ++q)
    {
      triSk1.push_back(*q);
    }
    triSk2.clear();
  }
}


void
adjustIndices(int locMap[3],TriSkeleton* triSk)
{
  int i,j;
  double pos[3][3];

  // set positions
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      pos[locMap[i]][j] = triSk->pos[i][j];
    }
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
       triSk->pos[i][j] = pos[i][j];
    }
  }

  // set children
  if (triSk->isRefined())
  {
    TriSkeleton* tmpCh[3];
    tmpCh[locMap[0]] =  triSk->child[0];
    tmpCh[locMap[1]] =  triSk->child[1];
    tmpCh[locMap[2]] =  triSk->child[2];

    triSk->child[0] = tmpCh[0];
    triSk->child[1] = tmpCh[1];
    triSk->child[2] = tmpCh[2];

    adjustIndices(locMap,triSk->child[0]);
    adjustIndices(locMap,triSk->child[1]);
    adjustIndices(locMap,triSk->child[2]);
    adjustIndices(locMap,triSk->child[3]);
  }
}


//
//
void
MultilevelGrid::setBoundaryMesh(vector<Element*>& element)
{
  list<Element*>* lt = grid[0]->getTriangleList();
  vector<float*> indexMap(lt->size());
  typedef list<Element*>::iterator TriI;
  for (TriI p = lt->begin(); p != lt->end(); ++p)
  {
    Element* tri = *p;
    Vertex* vt[3];
    vt[0] = tri->getVertex(0);
    vt[1] = tri->getVertex(1);
    vt[2] = tri->getVertex(2);
    int locMap[3];
    for (uint n = 0; n < element.size(); n++)
    {
      Vertex* chvt[3];
      chvt[0] = element[n]->getVertex(0);
      chvt[1] = element[n]->getVertex(1);
      chvt[2] = element[n]->getVertex(2);
      bool flag[3];
      flag[0] = false;
      flag[1] = false;
      flag[2] = false;
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          double dist = grd::distance(vt[i]->getPosition(),chvt[j]->getPosition());
          if (dist < 0.001)
          {
            locMap[i] = j;
            flag[i] = true;
            break;
          }
        }
      }
      if (flag[0] && flag[1] && flag[2]) {
        ((Triangle*)tri)->setSkeleton(((Triangle*)element[n])->getSkeleton());
        bool locFlag = false;
        for (int m = 0; m < 3; m++)
        {
          if (locMap[m] != m)
            locFlag = true;
        }
        if (locFlag == true) {
          adjustIndices(locMap,((Triangle*)element[n])->getSkeleton());
        }
        break;
      }
    }
  } // for loop over triangles
}


// Description
// top-down / bottom-up algorithm
//void
//MultilevelGrid::refine(int level,Stack<vector<int>*>& refElements,Stack<vector<int>*>& irregElements)
//{
//  int i;
//  int maxNoOfLevels;
//  int refLevel;
//  GridLevel* gridLevel;
//
//  int counter;
//  Cell* theCell;
//  CellListIterator iter;
//  Edge*    theEdge;
//  Element* theElement;
//
//  if (level < topLevel) {
//    for (i = topLevel; i > level; i--) {
//      gridLevel = grid.pop();
//      delete gridLevel;
//      noOfVertices.pop();
//      noOfEdges.pop();
//      noOfElements.pop();
//    }
//    topLevel = level;
//    grid[topLevel]->resetDescendant();
//  }
//  else if (level == topLevel)
//    return ;
//  else if (level > topLevel) {
//    maxNoOfLevels = refElements.size();
//    refLevel = (level < maxNoOfLevels) ? level : maxNoOfLevels;
//
//    for (i = topLevel; i < refLevel; i++) {
//      GridLevel* newGridLevel = new GridLevel;
//      // update data
//      grid.push(newGridLevel);
//      noOfVertices.push(noOfVertices[i]);
//      noOfEdges.push(noOfEdges[i]);
//      noOfElements.push(noOfElements[i]);
//
//      grid[i]->refine(grid[i+1],*refElements[i],*irregElements[i]);
//    }
//
//    topLevel = refLevel;
//    updateBoundary();
//
//  }
//
//  // recompute number of vertices
//  for (i = 1; i <= topLevel; i++)
//    noOfVertices[i]    = noOfVertices[i-1] + grid[i]->noOfVertices();
//
//  // recompute no. of edges
//  counter = 0;
//  noOfEdges[0] = grid[0]->noOfEdges();
//  for (i = 0; i < topLevel; i++) {
//    iter.init(grid[i]->getEdgeList());
//    while ((theCell = ++iter)) {
//      theEdge = (Edge*) theCell;
//      if (!(theEdge->isRefined()))
//        counter++;
//    }
//    noOfEdges[i+1] = counter + grid[i+1]->noOfEdges();
//  }
//
//  // recompute no. of elements
//  counter = 0;
//  noOfElements[0] = grid[0]->noOfElements();
//  for (i = 0; i < topLevel; i++) {
//    iter.init(grid[i]->getElementList());
//    while ((theCell = ++iter)) {
//      theElement = (Element*) theCell;
//      if (!(theElement->isRefined()))
//        counter++;
//    }
//    noOfElements[i+1] = counter + grid[i+1]->noOfElements();
//  }
//
//  // recompute no. of bnd. elements
//  counter = 0;
//  noOfBndElements[0] = grid[0]->noOfBndElements();
//  for (i = 0; i < topLevel; i++) {
//    iter.init(grid[i]->getBoundary());
//    while ((theCell = ++iter)) {
//      theElement = (Element*) theCell;
//      if (!(theElement->isRefined()))
//        counter++;
//    }
//    noOfBndElements[i+1] = counter + grid[i+1]->noOfBndElements();
//  }
//
//  bndTopLevel = topLevel;
//}






// Description
// Write compressed information for reconstructing
// the topology of the mesh hierarchy.
// This function does not write down the coarse mesh
//void
//MultilevelGrid::writeMLGrid(const char* filename,
//                              Vector<List<Cell*>*>& elemList,
//                              Vector<List<Cell*>*>& vertices,
//                              Vector<int>& vertexId)
//{
//  int i,j;
//  int ofd;
//  //int length;
//  int counter;
//
//  Vertex*  theVertex;
//  Edge*    theEdge;
//  Element* theElement;
//  Element* theChildElem;
//
//  CellListIterator      eIter;
//  ListIterator<Cell*> iter;
//  Cell**                theCellPtr;
//
//  // create vertex list
//  vertices.init(topLevel+1);
//
//  // create element lists
//  elemList.init(topLevel+1);
//
//  for (i = 0; i <= topLevel; i++) {
//    vertices[i] = new List<Cell*>;
//    elemList[i] = new List<Cell*>;
//  }
//
//
//  // first sort mesh
//  // used edge based refinement information
//  grid[0]->copyElementList(*elemList[0]);
//  grid[0]->copyVerticesList(*vertices[0]);
//
//  for (i = 0; i < topLevel; i++) {
//    iter.init(*elemList[i]);
//    while((theCellPtr = ++iter)) {
//      theElement = (Element*)(*theCellPtr);
//      if (theElement->isRefined()) {
//        // *- tetras
//         if ((theElement->getElementType()) == TETRAHEDRON) {
//           // first process edges
//           for (j = 0; j < 6; j++) {
//             theEdge = theElement->getEdge(j);
//             if (theEdge->isNotUsed()) {
//               theVertex = theEdge->getEdgeMidPoint();
//               vertices[i+1]->append(theVertex);
//               theEdge->setUsed();
//             }
//           }
//           // process tetra children
//           for (j = 0; j < 4; j++) {
//             theChildElem = theElement->getChild(j);
//             elemList[i+1]->append(theChildElem);
//           }
//           // process octa child, don't forget barycenter
//           theChildElem = theElement->getChild(4);
//           theVertex = theChildElem->getVertex(6);
//           vertices[i+1]->append(theVertex);
//           elemList[i+1]->append(theChildElem);
//         }
//        // *- octas
//         else if ((theElement->getElementType()) == OCTAHEDRON) {
//           // first process edges
//           for (j = 0; j < 12; j++) {
//             theEdge = theElement->getEdge(j);
//             if (theEdge->isNotUsed()) {
//               theVertex = theEdge->getEdgeMidPoint();
//               vertices[i+1]->append(theVertex);
//               theEdge->setUsed();
//             }
//           }
//           // process octa children
//           for (j = 0; j < 6; j++) {
//             theChildElem = theElement->getChild(j);
//             theVertex = theChildElem->getVertex(6);
//             vertices[i+1]->append(theVertex);
//             elemList[i+1]->append(theChildElem);
//           }
//           // process tetra children
//           for (j = 6; j < 14; j++) {
//             theChildElem = theElement->getChild(j);
//             elemList[i+1]->append(theChildElem);
//           }
//         }
//      }
//    }
//  }
//
//  // unset edges
//  unsetEdges();
//
//  // sort vertices
//  //length = getNoOfVertices();
//  //vertexId.init(length);
//  counter = 0;
//  for (i = 0; i <= topLevel; i++) {
//    iter.init(*vertices[i]);
//    while((theCellPtr = ++iter)) {
//      theVertex = (Vertex*)(*theCellPtr);
//      vertexId[counter++] = theVertex->getId();
//      theVertex->setId(counter);
//    }
//  }
//
//  // *-----------------------
//  // write hierarchical mesh
//  // *-----------------------
//  if ((ofd = open(filename,O_WRONLY|O_CREAT,0664)) == -1) {
//    WarningMacro("can't open output file");
//    exit(1);
//  }
//
//  // write no. of levels
//  write(ofd,&topLevel,sizeof(int));
//
//  // create and write refinement pattern
//  for (i = 0; i < topLevel; i++) {
//    int counter = 0;
//    Stack<int> refElem;
//    Stack<int> irregElem;
//    iter.init(*elemList[i]);
//    while((theCellPtr = ++iter)) {
//      theElement = (Element*)(*theCellPtr);
//      if (theElement->isRefined())
//         refElem.push(counter);
//      else if (theElement->isIrregular())
//         irregElem.push(counter);
//      counter++;
//    }
//
//    // write level i
//    write(ofd,&i,sizeof(int));
//    // write refined elems
//    int  noOfElem;
//    int* values;
//    noOfElem = refElem.size();
//    values   = refElem.getDattr();
//    write(ofd,&noOfElem,sizeof(int));
//    write(ofd,values,noOfElem*sizeof(int));
//    // write irreg elems
//    noOfElem = irregElem.size();
//    values   = irregElem.getDattr();
//    write(ofd,&noOfElem,sizeof(int));
//    write(ofd,values,noOfElem*sizeof(int));
//  }
//
//}
//
//
//// Description
//// This function read the information necessary
//// to reconstruct a mesh hierarchy starting with
//// a coarse mesh and a set of integer records.
//// The coarse mesh is read from a file
//// and the mesh hierarchy from a different one.
//void
//MultilevelGrid::readMLGrid(const char* mesh0,const char* mesh1)
//{
//  int i,j;
//  int ifd;
//  int noOfVertices;
//  int noOfEdges;
//  int noOfElements;
//  int noOfBndElems;
//  int newTopLevel;
//  int clevel;
//  int counter;
//  int buff[20];
//
//  double pos[3];
//
//  CellListIterator iter;
//  BoundaryIterator bIter;
//  Element* theElement;
//
//  Vector<Vertex*>  vertex;
//  Vector<Element*> element;
//
//  Stack<Vector<int>*> refElements;
//  Stack<Vector<int>*> irregElements;
//
//  if ((ifd = open(mesh0,O_RDONLY)) == -1) {
//    WarningMacro("can't open input file");
//    exit(1);
//  }
//
//  // read no. of vertices
//  read(ifd,&noOfVertices,sizeof(int));
//  // read no. of edges
//  read(ifd,&noOfEdges,sizeof(int));
//  // read no. of elements
//  read(ifd,&noOfElements,sizeof(int));
//
//  // read no. of bnd. elements
//  read(ifd,&noOfBndElems,sizeof(int));
//
//  // read no. of face sets
//  read(ifd,&noOfFaceSets,sizeof(int));
//
//  // set vertex array
//  vertex.init(noOfVertices);
//  for (i = 0; i < noOfVertices; i++) {
//    Vertex* theVertex = new Vertex;
//    theVertex->setId(i+1);
//    vertex[i] = theVertex;
//  }
//
//  // read element to vertex connectivity
//  element.init(noOfElements);
//  for (i = 0; i < noOfElements; i++) {
//    read(ifd,buff,sizeof(int));
//    // *- tets
//    if (buff[0] == 4) {
//      Tetrahedron* tet = new Tetrahedron;
//      element[i] = tet;
//    }
//    read(ifd,&(buff[1]),buff[0]*sizeof(int));
//    for (j = 0; j < buff[0]; j++)
//      element[i]->setVertex(j,vertex[buff[j+1] - 1]);
//  }
//
//  // compute mesh topology
//  buildCoarseMesh(vertex,element);
//
//  // read vertex positions
//  for (i = 0; i < noOfVertices; i++) {
//    read(ifd,pos,3*sizeof(double));
//    vertex[i]->setPosition(pos);
//  }
//
//  // compute boundary
//  extractBoundary();
//
//  if (noOfBndElems != noOfBndElements[0]) {
//    ErrorMacro("wrong number of boundary elements");
//    exit(1);
//  }
//
//  // set element set tag
//  char* setTags = new char[noOfBndElems];
//  read(ifd,setTags,noOfBndElems*sizeof(int));
//
//  counter = 0;
//  bIter.init(*this);
//  while((theElement = ++ bIter))
//    theElement->setElementSet(setTags[counter++]);
//
//
//  // count vertices, edges and elements
//  initElemCount();
//
//  // contruct abstract description of multilevel mesh
//  if ((ifd = open(mesh1,O_RDONLY)) == -1) {
//    WarningMacro("can't open input file");
//    exit(1);
//  }
//
//  // read no. of levels
//  read(ifd,&newTopLevel,sizeof(int));
//
//  for (i = 0; i < newTopLevel; i++) {
//    read(ifd,&clevel,sizeof(int));
//    if (clevel != i) {
//      WarningMacro("reading wrong level");
//      exit(1);
//    }
//    int noOfRefElems;
//    int noOfIrregElems;
//    read(ifd,&noOfRefElems,sizeof(int));
//
//    Vector<int>* refElem = new Vector<int>(noOfRefElems);
//    read(ifd,refElem->getDattr(),noOfRefElems*sizeof(int));
//
//    read(ifd,&noOfIrregElems,sizeof(int));
//    Vector<int>* irregElem = new Vector<int>(noOfIrregElems);
//
//    read(ifd,irregElem->getDattr(),noOfIrregElems*sizeof(int));
//
//    refElements.push(refElem);
//    irregElements.push(irregElem);
//  }
//
//  refine(newTopLevel,refElements,irregElements);
//
//  setVertexNumbers();
//
//}
//
//



} // namespace grd
