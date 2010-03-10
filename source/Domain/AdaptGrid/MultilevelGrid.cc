// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/***************************************************************************
    File        : MultilevelGrid.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Jan 18 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "MultilevelGrid.hh"

// Utilities
#include "Utility.hh"

using namespace std;

namespace grd {


// Description
//
MultilevelGrid::MultilevelGrid()
{
  // Marks and tags
  closure  = true;
  boundary = true;
  // create grid level 0
  topLevel = 0;
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
  for (size_t i = 0; i < grid.size(); i++)
    delete grid[i];
}



// Description
// top-down / bottom-up algorithm
void
MultilevelGrid::refine()
{
  // Loop index
  int i;

  // Create new grid level and push into stack
  GridLevel* newGridLevel = new GridLevel;
  grid.push(newGridLevel);

  // Unsplit regular octahedra
  for (int level = 0; level <= topLevel; level++) grid[level]->unsplitOctahedra();
  for (int level = 0; level <= topLevel; level++) grid[level]->removeOctahedronSplitElements();

  // Check Element Marks for consistency
  evaluateRefinementMarks();

  // Prepare Mesh for Refinement
  prepareGridForRefinement();
  
  // top-down
  // refine marked elements
  // irregular elements are refined regular
  for (i = topLevel; i >= 0; i--) grid[i]->refine(*grid[i+1]);

  // bottom-up
  // coarsen grid level
  for(i = 0; i <= topLevel; i++) grid[i]->coarsen();

  // remove unused vertices and edges from lists
  //   1. reference vertices and edges
  //   2. remove unreferenced vertices and edges
  //   3. unreference vertices and edges
  for (i = 0; i <= (topLevel+1); i++) grid[i]->refVerticesAndEdges();
  for (i = 0; i <= (topLevel+1); i++) grid[i]->removeUnrefVerticesAndEdges();
  for (i = 0; i <= (topLevel+1); i++) grid[i]->unrefVerticesAndEdges();

  //close grid level
  // needs to number vertices
  setVertexNumbers();
  for(i = 0; i <= topLevel; i++) grid[i]->close(*grid[i+1]);

  // compute new no. of levels
  if (!grid[topLevel+1]->empty())
  {
    topLevel++;
  }
  else if (topLevel > 0)
  {
    delete grid[topLevel+1];
    grid.pop();
    if (grid[topLevel]->empty()) {
      delete grid[topLevel];
      grid.pop();
      topLevel--;
    }
  }

  // Complete Mesh
  if (closure) completeGridWithClosureElements();

  // Split octahedra
  for (int level = 0; level <= topLevel; level++) grid[level]->splitOctahedra();

  // Recompute number of elements
  updateNoOfElements();
}


// Description
//  refine all elements uniformly
void
MultilevelGrid::uniformRefine()
{
  // Mark all elements for refinement
  for (int level = 0; level <= topLevel; level++)
  {
    typedef std::list<Element*>::iterator ElmI;

    std::list<Element*>** lt = grid[level]->getElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        // call functor
        if (!(*p)->isRefined())
          (*p)->markForRefinement();
      }
      // next element type
      type++;
    }
  }

  // Refine Grid
  refine();
}



//
// Description
// compute mesh topology from a list of vertices and a list of elements
void
MultilevelGrid::buildCoarseGrid(vector<Vertex*>& vertex,vector<Element*>& element)
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
  list<Element*>* hexList    = grid[0]->getHexahedronList();
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
      case GRD_HEXAHEDRON:
        hexList->push_back(element[i]);
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
      int noOfVerticesAtFace = theElement->noOfVerticesAtFace(face);
      Stack<Element*>** stck = new Stack<Element*>*[noOfVerticesAtFace];
      for (m = 0; m < noOfVerticesAtFace; m++)
      {
        vt = theElement->getVertexAtFace(face,m);
        stck[m] = vtList[vt->getId() - 1];
      }

      // Init:
      // intersection of first two lists
      list<Element*>* lst1 = new list<Element*>;
      list<Element*>* lst2 = new list<Element*>;
      for (i = 0; i < (int)stck[0]->size(); i++)
      {
        neigh1 = (*stck[0])[i];
        if (neigh1 != theElement)
        {
          for (j = 0; j < (int)stck[1]->size(); j++)
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
      for (m = 2; m < noOfVerticesAtFace; m++)
      {
        // intersect allways two lists
        for (i = 0; i < (int)stck[m]->size(); i++)
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
        WARN("more than one neighbor for a face");

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
    if (setCounter > 1)
    {
      Error("panic! set counter > 1");
    }
    for (edge = 0; edge < theElement->getNoOfEdges(); edge++) {
      setCounter = 0;
      theEdge = theElement->getEdge(edge);
      if (!theEdge) {
        vt1 = theElement->getVertexAtEdge(edge,0);
        stack1 = vtList[vt1->getId() - 1];
        vt2 = theElement->getVertexAtEdge(edge,1);
        stack2 = vtList[vt2->getId() - 1];

        for (i = 0; i < (int)stack1->size(); i++) {
          if (!(theElement->getEdge(edge))) {
            neigh1 = (*stack1)[i];
            if (neigh1 != theElement) {
              for (j = 0; j < (int)stack2->size(); j++) {
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
  if (isVolume() && boundary)
    extractBoundary();

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
  int size;
  int noOfLevels = topLevel+1;

  // First have to update stack sizes
  size = noOfVertices.size();
  if (size < noOfLevels)
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
  }
  else if (size > noOfLevels)
  {
    noOfVertices.pop();
    noOfEdges.pop();
    noOfTriangles.pop();
    noOfQuadrangles.pop();
    noOfTetrahedra.pop();
    noOfOctahedra.pop();
    noOfHexahedra.pop();
    noOfPrisms.pop();
    noOfPyramids.pop();
  }
  
  // VERTICES
  // Update the no. of vertices for each level
  size = noOfVertices.size();
  noOfVertices[0] = grid[0]->noOfVertices();
  for (i = 1; i <= topLevel; i++)
    noOfVertices[i]    = noOfVertices[i-1] + grid[i]->noOfVertices();

  // EDGES
  // Update the no. of edges for each level
  typedef list<Edge*>::iterator EdgI;
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfEdges[i] = grid[i]->noOfEdges() + counter;                    
    list<Edge*>* lt = grid[i]->getEdgeList();
    for (EdgI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!((*p)->isRefined()))
        counter++;
    }
  }
  noOfEdges[topLevel] = counter + grid[topLevel]->noOfEdges();

  // Compute no. of Elements
  // Iterate over lists of elements
  typedef list<Element*>::iterator EI;
  list<Element*>* lt;

  // TRIANLGES
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfTriangles[i] = grid[i]->noOfTriangles() + counter;
    // Triangles
    lt = grid[i]->getTriangleList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if ((*p)->isRegular() || (*p)->isClosure())
        counter++;
    }
  }
  noOfTriangles[topLevel] = counter + grid[topLevel]->noOfTriangles();
  // QUADRANGLES
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfQuadrangles[i] = grid[i]->noOfQuadrangles();
    // Triangles
    lt = grid[i]->getQuadrangleList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
  }
  noOfQuadrangles[topLevel] = counter + grid[topLevel]->noOfQuadrangles();
  // Tetrahedra
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfTetrahedra[i] = grid[i]->noOfTetrahedra();
    // Triangles
    lt = grid[i]->getTetrahedronList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
  }
  noOfTetrahedra[topLevel] = counter + grid[topLevel]->noOfTetrahedra();
  // Octahedra
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfOctahedra[i] = grid[i]->noOfOctahedra();
    // Triangles
    lt = grid[i]->getOctahedronList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
  }
  noOfOctahedra[topLevel] = counter + grid[topLevel]->noOfOctahedra();
  // Hexahedra
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfHexahedra[i] = grid[i]->noOfHexahedra();
    // Triangles
    lt = grid[i]->getHexahedronList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
  }
  noOfHexahedra[topLevel] = counter + grid[topLevel]->noOfHexahedra();
  // Prisms
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfPrisms[i] = grid[i]->noOfPrisms();
    // Triangles
    lt = grid[i]->getPrismList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
  }
  noOfPrisms[topLevel] = counter + grid[topLevel]->noOfPrisms();
  // Pyramids
  counter = 0;
  for (i = 0; i < topLevel; i++) 
  {
    noOfPyramids[i] = grid[i]->noOfPyramids();
    // Triangles
    lt = grid[i]->getPyramidList();
    for (EI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        counter++;
    }
  }
  noOfPyramids[topLevel] = counter + grid[topLevel]->noOfPyramids();
}


// Description
//  Compute the boundary of a volumen mesh
void
MultilevelGrid::extractBoundary()
{
  int i;
  int face;
  Element* theElement;
  Element* theNeighbor;

  // Create lists of Triangles and Quadrilaterals
  std::list<Element*>* triList = getTriangleList(0);
  std::list<Element*>* quaList = getQuadrangleList(0);

  typedef std::list<Element*>::iterator ElmI;
  std::list<Element*>** lt = grid[0]->getElementList();
  int type = 2;
  while (lt[type])  // for each element type
  {
    for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
    {
      theElement = *p;
      for (face = 0; face < theElement->getNoOfFaces(); face++)
      {
        theNeighbor = theElement->getNeighbor(face);
        if (!theNeighbor)
        {
          int noOfVertex = theElement->noOfVerticesAtFace(face);
          // The neighbor is a triangle
          if (noOfVertex == 3)
          {
            Triangle* tri = new Triangle;
            Vertex* vt[3];
            Edge*   eg[3];
            for (i = 0; i < 3; i++)
            {
              vt[i] = theElement->getVertexAtFace(face,i);
              eg[i] = theElement->getEdgeAtFace(face,i);
            }
            // Assemble triangle
            Vertex* vv1;
            Vertex* vv2;
            // set edges
            for (int edge = 0; edge < 3; edge++)
            {
              tri->setEdge(edge,eg[edge]);
              vv1 = eg[edge]->getVertex(0);
              vv2 = eg[edge]->getVertex(1);
              for (i = 0; i < 3; i++)
              {
                if (vt[i] != vv1 && vt[i] != vv2)
                  tri->setVertex(edge,vt[i]);
              }
            }
            triList->push_back(tri);
            tri->setNeighbor(0,theElement);
            theElement->setNeighbor(face,tri);
          }
          // The neighbor is a quadrangle
          else if (noOfVertex == 4)
          {
            Quadrangle* qua = new Quadrangle;
            Vertex* vt[4];
            Edge*   eg[4];
            for (i = 0; i < 4; i++)
            {
              vt[i] = theElement->getVertexAtFace(face,i);
              eg[i] = theElement->getEdgeAtFace(face,i);
            }
            // Set quad
            Vertex* v0;
            Vertex* v1;
            // edge 0
            v0 = eg[0]->getVertex(0);
            v1 = eg[0]->getVertex(1);
            if ((v0 == vt[0] && v1 == vt[1]) || (v0 == vt[1] && v1 == vt[0]))
            {
              qua->setVertex(0,vt[0]);
              qua->setEdge(0,eg[0]);
            } else std::cout << "ERROR: wrong edge 0" << std::endl;
            v0 = eg[1]->getVertex(0);
            v1 = eg[1]->getVertex(1);
            if ((v0 == vt[1] && v1 == vt[2]) || (v0 == vt[2] && v1 == vt[1]))
            {
              qua->setVertex(1,vt[1]);
              qua->setEdge(1,eg[1]);
            } else std::cout << "ERROR: wrong edge 1" << std::endl;
            v0 = eg[2]->getVertex(0);
            v1 = eg[2]->getVertex(1);
            if ((v0 == vt[2] && v1 == vt[3]) || (v0 == vt[3] && v1 == vt[2]))
            {
              qua->setVertex(2,vt[2]);
              qua->setEdge(2,eg[2]);
            } else std::cout << "ERROR: wrong edge 2" << std::endl;
            v0 = eg[3]->getVertex(0);
            v1 = eg[3]->getVertex(1);
            if ((v0 == vt[3] && v1 == vt[0]) || (v0 == vt[0] && v1 == vt[3]))
            {
              qua->setVertex(3,vt[3]);
              qua->setEdge(3,eg[3]);
            } else std::cout << "ERROR: wrong edge 3" << std::endl;
            
            quaList->push_back(qua);
            qua->setNeighbor(0,theElement);
            theElement->setNeighbor(face,qua);
          }
        } // the neighbor is null
      } // for each element's face
    } // for each element

    // get next element type
    type++;
  } // while element type
}


void
MultilevelGrid::evaluateRefinementMarks()
{
  // Elements of the coarse mesh cannot be coarsened
  typedef std::list<Element*>::iterator ElmI;
  std::list<Element*>** lt = grid[0]->getElementList();
  int type = 0;
  while (lt[type])
  {
    for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
    {
      if ((*p)->isRegular() && (*p)->isMarkedForCoarsening())
      {
        (*p)->unmarkForCoarsening();
      }
    }
    // next element type
    type++;
  }

  // Check refinement marks for boundary elements
  for (int j = 0; j <= topLevel; j++)
  {
    std::list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type] && type < 2)
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRefined())
          continue;
        if (((*p)->type() == GRD_TRIANGLE) || ((*p)->type() == GRD_QUADRANGLE))
        {
          // the neighbor volume element
          Element* neighbor = (*p)->getNeighbor(0);
          if (neighbor)
          {
            if (neighbor->isMarkedForRefinement())
              (*p)->markForRefinement();
            else if (neighbor->isMarkedForCoarsening())
              (*p)->markForCoarsening();
            else if ((*p)->isMarkedForRefinement())
              neighbor->markForRefinement();
            else if ((*p)->isMarkedForCoarsening())
              neighbor->markForCoarsening();
          }
        }
      }
      // next element type
      type++;
    }
  }

  // Mark irregular elements for ref
  // or coarsening
  // Refined elements cannot be marked for refinement
  for (int j = 0; j <= topLevel; j++)
  {
    std::list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRefined() && ((*p)->isMarkedForRefinement() || (*p)->isMarkedForCoarsening()))
        {
          (*p)->unsetAllRefinementMarks();
        }
        else if ((*p)->isClosure())
        {
          // the neighbor volume element
          if ((*p)->isMarkedForRefinement())
          {
            (*p)->markParentForRefinement();
          }
        }
      }
      // next element type
      type++;
    }
  }
}


/////////////////////////////////////////////////////////////////////
// Auxiliary Functions
/////////////////////////////////////////////////////////////////////
inline bool compareVerticesAtFace(Vertex* vt[3],Vertex* vv[3],int nov);
inline bool checkElementConnectivity(Element* t);
inline bool checkElementTopology(Element* t);

// Description
//  Remove closure elements
//  The refine algorithm works only for meshes
//  with irregular elements
void
MultilevelGrid::prepareGridForRefinement()
{
  typedef std::list<Element*>::iterator ElmI;
  // Check refinement marks for boundary elements
  for (int level = 0; level <= topLevel; level++)
  {
    std::list<Element*>** lt = grid[level]->getElementList();
    int type = 0;
    while (lt[type])
    {
      ElmI p = lt[type]->begin();
      while (p != lt[type]->end())
      {
          if ((*p)->isIrregular())
          {
            for (int nn = 0; nn < (*p)->getNoOfFaces(); nn++)
            {
              Element* neighbor = (*p)->getNeighbor(nn);
              if (neighbor)
              {
                if (neighbor->isClosure())
                {
                  (*p)->setNeighbor(nn,0);
                }
                else
                {
                  Vertex* vt[4];
                  Vertex* vv[4];
                  int noOfVt = (*p)->noOfVerticesAtFace(nn);
                  for (int i = 0; i < noOfVt; i++) vt[i] = (*p)->getVertexAtFace(nn,i);
                  for (int fc = 0; fc < neighbor->getNoOfFaces(); fc++)
                  {
                    int noOfVv = neighbor->noOfVerticesAtFace(fc);
                    if (noOfVt == noOfVv)
                    {
                      for (int j = 0; j < noOfVt; j++) vv[j] = neighbor->getVertexAtFace(fc,j);
                      if (compareVerticesAtFace(vt,vv,noOfVt))
                      {
                        // neighbour found
                        neighbor->setNeighbor(fc,*p);
                        break;
                      }
                    }
                  }
                } // else neighbor is not a triangle
              } // if (neighbor)
            } // for noOfFaces of p

            // reset the pointer to the children
            (*p)->resetPointerToChildren();
          }
          // update pointer for the while loop
          ++p;
      } // while p != lt[type]->end()
      // next element type
      type++;
    }
  }



  for (int level = 0; level <= topLevel; level++)
  {
    std::list<Element*>** lt = grid[level]->getElementList();
    int type = 0;
    while (lt[type])
    {
      ElmI p = lt[type]->begin();
      while (p != lt[type]->end())
      {
        if ((*p)->isClosure())
        {
          for (int i = 0; i < (*p)->getNoOfFaces(); i++)
          {
            Element* neighbor = (*p)->getNeighbor(i);
            if (neighbor && !neighbor->isClosure())
            {
              int cface = neighbor->getCommonFace(*p);
              if (cface != -1)
                neighbor->setNeighbor(cface,0);
            }
          }
          // delete element pointed to
          delete *p;
          // remove pointer from the list
          p = lt[type]->erase(p);
        } // if p is closure
        else
          ++p;
      } // while p != lt[type]->end()
      // next element type
      type++;
    }
  }


  // Remove unused elements
  // remove unused vertices and edges from lists
  //   1. reference vertices and edges
  //   2. remove unreferenced vertices and edges
  //   3. unreference vertices and edges
  for (int i = 0; i <= topLevel; i++) grid[i]->refVerticesAndEdges();
  for (int i = 0; i <= topLevel; i++) grid[i]->removeUnrefVerticesAndEdges();
  for (int i = 0; i <= topLevel; i++) grid[i]->unrefVerticesAndEdges();

}

//
//
// Compute the conforming clousre
void
MultilevelGrid::completeGridWithClosureElements()
{
  typedef std::list<Element*>::iterator ElmI;
  for (int j = 0; j < topLevel; j++)
  {
    std::list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isIrregular())
        {
          (*p)->close(*grid[j+1]);
        }
      }
      // next element type
      type++;
    }
  }
}




// Controll mesh topology and connectivity
bool
MultilevelGrid::checkGridConsistencyAndConnectivity()
{
  typedef std::list<Element*>::iterator ElmI;
  for (int j = 0; j <= topLevel; j++)
  {
    std::list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type])
    {
      if (type == 3)
      {
        // check neighbors before unsplit
        for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
        {
          if ((*p)->isRegular())
          {
            for (int nn = 0; nn < 4; nn++)
            {
              Element* child = (*p)->getChild(nn);
              for (int mm = 0; mm < 4; mm++)
              {
                Element* neighbor = child->getNeighbor(mm);
                int cface = neighbor->getCommonFace(child);
                Element* myself = neighbor->getNeighbor(cface);
                if (myself != child)
                  std::cout << "ERROR: wrong connectivity" << std::endl;
              }
            }
          }
        }         
        type++;
        continue;
      }
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRegular() || (*p)->isClosure())
        {
          if (!checkElementConnectivity(*p))
          {
            std::cout << "ERROR: wrong element connectivity" << std::endl;
            return false;
          }
        }
        // check edges here
        for (int nn = 0; nn < (*p)->getNoOfEdges(); nn++)
        {
          Edge* ee = (*p)->getEdge(nn);
          if (!ee)
          {
            std::cout << "ERROR: element has a zero edge, type " << ((*p)->type()) << std::endl;
          }
        }
      }
      // next element type
      type++;
    }
  }

  // finish
  return true;
}


//
// implement the auxiliary functions
//
inline bool compareVerticesAtFace(Vertex* vt[3],Vertex* vv[3],int nov)
{
  bool flag = false;
  int count = 0;
  for (int i = 0; i < nov; i++)
  {
    for (int j = 0; j < nov; j++)
    {
      if (vt[i] == vv[j])
      {
        count++;
        break;
      }
    }
  }
  if (count == nov)
    flag = true;

  return flag;
}

// check elements for connectivity
inline bool checkElementConnectivity(Element* t)
{
  Vertex* vt[3];
  Vertex* vv[3];
  // Sort elements by type
  if (t->type() == GRD_TRIANGLE)
  {
    // check connection with the volume
    // element has only one "face" and three vertices
    vt[0] = t->getVertex(0);
    vt[1] = t->getVertex(2);
    vt[2] = t->getVertex(1);

    // get neighbor
    Element* neighbor = t->getNeighbor(0);
    if (!neighbor)
    {
      std::cout << "ERROR: triangle has no neighbor" << std::endl;
      return false;
    }
    for (int nn = 0; nn < neighbor->getNoOfFaces(); nn++)
    {
      vv[0] = neighbor->getVertexAtFace(nn,0);
      vv[1] = neighbor->getVertexAtFace(nn,1);
      vv[2] = neighbor->getVertexAtFace(nn,2);
      if (compareVerticesAtFace(vt,vv,3))
      {
        return true;
      }
    }
  }
  else if (t->type() == GRD_TETRAHEDRON || t->type() == GRD_OCTAHEDRON)
  {
    int counter = 0;
    for (int nn = 0; nn < t->getNoOfFaces(); nn++)
    {
      vt[0] = t->getVertexAtFace(nn,0);
      vt[1] = t->getVertexAtFace(nn,1);
      vt[2] = t->getVertexAtFace(nn,2);
      Element* neighbor = t->getNeighbor(nn);
      if (neighbor)
      {
        if (neighbor->type() == GRD_TRIANGLE)
        {
          vv[0] = neighbor->getVertex(0);
          vv[1] = neighbor->getVertex(1);
          vv[2] = neighbor->getVertex(2);
          if (compareVerticesAtFace(vt,vv,3))
          {
            counter++;
          }
        }
        else
        {
          for (int mm = 0; mm < neighbor->getNoOfFaces(); mm++)
          {
            vv[0] = neighbor->getVertexAtFace(mm,0);
            vv[1] = neighbor->getVertexAtFace(mm,1);
            vv[2] = neighbor->getVertexAtFace(mm,2);
            if (compareVerticesAtFace(vt,vv,3))
            {
              counter++;
              break;
            }
          }
        }

      }
      else
      {
        // catastrophic error
        std::cout << "ERROR: element has no neighbor" << std::endl;
        return false;
      }
    }
    // check for all neighbors
    if (counter == t->getNoOfFaces())
      return true;
    else
    {
      std::cout << "ERROR: can't find common face with neighbor" << std::endl;
      return false;
    }
  }

  // Compare pointers to neighbors
  for (int nn = 0; nn < t->getNoOfFaces(); nn++)
  {
    Element* neighbor = t->getNeighbor(nn);
    int cface = neighbor->getCommonFace(t);
    Element* myself = neighbor->getNeighbor(cface);
    if (myself != t)
      std::cout << "ERROR: wrong pointer to neighbors" << std::endl;
  }
  
  // finish
  return false;
}


inline bool checkElementTopology(Element* t)
{
  for (int nn = 0; nn < t->getNoOfEdges(); nn++)
  {
    Edge* ee = t->getEdge(nn);
    if (!ee)
    {
      return false;
    }
  }
  return true;
}



//
void
MultilevelGrid::checkEdgesOfAllElements(char* text)
{
  typedef std::list<Element*>::iterator ElmI;
  // check all edges now
  for (int level = 0; level <= topLevel; level++)
  {
    std::list<Element*>** lt = grid[level]->getElementList();
    int type = 0;
    while (lt[type])
    {
      for  (ElmI q = lt[type]->begin(); q != lt[type]->end(); ++q)
      {
        if (checkElementTopology(*q) == false)
          std::cout << "ERROR: found for element: " << ((*q)->type()) << "  at " << text << std::endl;
      }
      type++;
    } 
  }
  
  // check all edges
  for (int level = 0; level <= topLevel; level++)
  {
    std::list<Edge*>* lt = grid[level]->getEdgeList();
    std::list<Edge*>::iterator q;
    for (q = lt->begin(); q != lt->end(); ++q)
    {
      if ((*q)->queryRef())
        std::cout << "ERROR: Edge is referenced" << std::endl;
    }
  }
}



//===========================================================================================
// Data Access
//===========================================================================================
//===========================================================================================

void
MultilevelGrid::getVolumeGrid(std::vector<Element*>& gr)
{
  gr.clear();
  typedef std::list<Element*>::iterator ElmI;
  for (int jj = 0; jj <= topLevel; jj++)
  {
    std::list<Element*>** lt = grid[jj]->getElementList();
    int type = 2;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if (!(*p)->isRefined())
        {
          gr.push_back(*p);
        }
      }
      // next element type
      type++;
    }
  }
}

// // Check quads and hexs
// void
// MultilevelGrid::checkQuadsAndHexs()
// {
//   for (int level = 0; level <= topLevel; level++)
//   {
//     std::list<Element*>** lt = grid[level]->getElementList();
//     std::list<Element*>::iterator q;
//     int typeQua = 1;
//     int typeHex = 4;
//     int faceNo = 0;
//     for  (q = lt[typeQua]->begin(); q != lt[typeQua]->end(); ++q)
//     {
//       std::cout << " ... processing face: " << faceNo << "  at level: " << level << std::endl;
//       Element* hex = (*q)->getNeighbor(0);
//       Vertex* qvt[4];
//       Vertex* hvt[8];
//       int noOfVert = hex->getNoOfVertices();
//       for (int nn = 0; nn < noOfVert; nn++) hvt[nn] = hex->getVertex(nn);
//       noOfVert = (*q)->getNoOfVertices();
//       for (int nn = 0; nn < noOfVert; nn++) qvt[nn] = (*q)->getVertex(nn);
//       int counter = 0;
//       for (int nn = 0; nn < 4; nn++)
//       {
//         for (int mm = 0; mm < 8; mm++)
//           if (qvt[nn] == hvt[mm]) { counter++; break;}
//       }
//       if (counter != 4) std::cout << "ERROR: quad's vertices not in neighbor hex with: " << counter <<  std::endl;
// 
//       Edge* qeg[4];
//       Edge* heg[12];
//       for (int nn = 0; nn < 12; nn++) heg[nn] = hex->getEdge(nn);
//       noOfVert = (*q)->getNoOfVertices();
//       for (int nn = 0; nn <  4; nn++) qeg[nn] = (*q)->getEdge(nn);
//       counter = 0;
//       for (int nn = 0; nn < 4; nn++)
//       {
//         for (int mm = 0; mm < 12; mm++)
//           if (qeg[nn] == heg[mm]) { counter++; break;}
//       }
//       if (counter != 4) std::cout << "ERROR: quad's edges not in neighbor hex with: " << counter <<  std::endl;
//       
//       // Check edges
//       int face = hex->getCommonFace(*q);
//       for (int nn = 0; nn < (*q)->getNoOfEdges(); nn++)
//       {
//         if (!((*q)->getEdge(nn))) std::cout << "ERROR: quad edge " << nn << " at hex face " << face << " is NULL" << std::endl;
//       }
//       faceNo++;
//     }
//     for  (q = lt[typeHex]->begin(); q != lt[typeHex]->end(); ++q)
//     {
//       // Check edges
//       for (int nn = 0; nn < (*q)->getNoOfEdges(); nn++)
//       {
//         if (!((*q)->getEdge(nn))) std::cout << "ERROR: hex edge " << nn << " is NULL" << std::endl;
//       }
//     }
//   }
// }

} // namespace grd
